#include "conforming_refinement.h"
#include "pillow_geometry.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <functional>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>

namespace refinement {
namespace {
constexpr int DEN=24*24*24;
constexpr int BITS[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                        {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
constexpr int FACES[6][4]={{0,3,2,1},{4,5,6,7},{0,1,5,4},
                          {1,2,6,5},{2,3,7,6},{3,0,4,7}};
// Face coordinates match compute_face_signature: xy, xz, or yz.
constexpr int BASIS[6][4]={{0,1,2,3},{4,5,6,7},{0,1,5,4},
                          {1,2,6,5},{3,2,6,7},{0,3,7,4}};
using Quad=std::array<int,4>;
using Weights=std::vector<std::pair<int,int>>;
using Tile=std::array<RationalCoord2,4>;
using Signature=std::vector<Tile>;
template<class T> std::array<T,4> ring(std::array<T,4> q) {
    auto best=q;
    for(int rev=0;rev<2;++rev) for(int shift=0;shift<4;++shift) {
        std::array<T,4> v;
        for(int k=0;k<4;++k) v[k]=q[(shift+(rev ? -k:k)+8)%4];
        best=std::min(best,v);
    }
    return best;
}
HexCorners geometry(const MeshData &m,int ie) {
    HexCorners h;
    for(int k=0;k<8;++k) h[k]=m.points.at(m.cells[ie][k]);
    return h;
}
std::array<RationalCoord3,8> unit_points() {
    std::array<RationalCoord3,8> p;
    for(int k=0;k<8;++k) p[k]={int8_t(24*BITS[k][0]),int8_t(24*BITS[k][1]),int8_t(24*BITS[k][2])};
    return p;
}
Weights weights(const HexConnectivity &cell,const RationalCoord3 &p) {
    std::map<int,int> w;
    int q[3]={p.x,p.y,p.z};
    for(int k=0;k<8;++k) {
        int a=1;
        for(int d=0;d<3;++d) a*=BITS[k][d] ? q[d]:24-q[d];
        if(a) w[cell[k]]+=a;
    }
    return {w.begin(),w.end()};
}
struct Option {
    int pad,children;
    std::vector<std::array<RationalCoord3,8>> cells;
    std::array<Signature,6> faces;
};
std::vector<Option> options() {
    const auto u=unit_points();
    HexCorners cube;
    for(int k=0;k<8;++k) cube[k]={u[k].x_double(),u[k].y_double(),u[k].z_double()};
    Option unchanged{0,1,{u},{}};
    for(int f=0;f<6;++f) unchanged.faces[f]={Tile{{{0,0},{24,0},{24,24},{0,24}}}};
    std::vector<Option> out={unchanged};
    for(int pad:get_all_pad_codes()) {
        Option o; o.pad=pad;
        auto children=generate_children_for_pad(pad,cube);
        o.children=children.size();
        for(const auto &h:children) {
            std::array<RationalCoord3,8> c;
            for(int k=0;k<8;++k) c[k]={int8_t(std::lround(h[k].x*24)),
                int8_t(std::lround(h[k].y*24)),int8_t(std::lround(h[k].z*24))};
            o.cells.push_back(c);
        }
        for(int f=0;f<6;++f) o.faces[f]=compute_face_signature_for_pad(pad,f).quads;
        out.push_back(std::move(o));
    }
    std::stable_sort(out.begin(),out.end(),[](const Option &a,const Option &b){
        return std::make_pair(a.children,a.pad)<std::make_pair(b.children,b.pad);
    });
    return out;
}
Signature canonical_face(const Signature &s,const HexConnectivity &cell,int f) {
    Quad ids;
    for(int k=0;k<4;++k) ids[k]=cell[BASIS[f][k]];
    const Quad canonical=ring(ids);
    const RationalCoord2 uv[4]={{0,0},{24,0},{24,24},{0,24}};
    RationalCoord2 target[4];
    for(int i=0;i<4;++i) for(int j=0;j<4;++j)
        if(ids[i]==canonical[j]) target[i]=uv[j];
    Signature out=s;
    for(auto &quad:out) {
        for(auto &p:quad) {
            const int x=p.u,y=p.v;
            p={int8_t(target[0].u+((target[1].u-target[0].u)*x+(target[3].u-target[0].u)*y)/24),
               int8_t(target[0].v+((target[1].v-target[0].v)*x+(target[3].v-target[0].v)*y)/24)};
        }
        quad=ring(quad);
    }
    std::sort(out.begin(),out.end());
    return out;
}
} // namespace

RefinedMesh refine_conforming(const MeshData &mesh,const std::vector<int> &requests,
                             std::size_t max_elements,std::size_t max_steps) {
    if(mesh.cells.size()>max_elements) throw std::runtime_error("Refinement exceeds refinementMaxElements before insertion");
    for(const auto &c:mesh.cells) {
        std::set<int> distinct(c.begin(),c.end());
        if(distinct.size()!=8) throw std::runtime_error("Refinement parent repeats a node");
        for(int id:c) if(id<0 || std::size_t(id)>=mesh.points.size()) throw std::runtime_error("Refinement parent node out of range");
    }
    static const auto opt=options();
    const int n=mesh.cells.size(), no=opt.size();
    std::map<Quad,std::vector<std::pair<int,int>>> owners;
    for(int i=0;i<n;++i) for(int f=0;f<6;++f) {
        Quad q;for(int k=0;k<4;++k)q[k]=mesh.cells[i][FACES[f][k]];
        std::sort(q.begin(),q.end());owners[q].push_back({i,f});
    }
    struct Link {int a,fa,b,fb;};
    std::vector<Link> links;
    for(auto &p:owners) {
        if(p.second.size()>2) throw std::runtime_error("Non-manifold parent face");
        if(p.second.size()==2) {
            auto a=p.second[0],b=p.second[1];
            Quad qa,qb;for(int k=0;k<4;++k){qa[k]=mesh.cells[a.first][FACES[a.second][k]];qb[k]=mesh.cells[b.first][FACES[b.second][k]];}
            if(ring(qa)!=ring(qb)) throw std::runtime_error("Parent faces have incompatible corner rings");
            links.push_back({a.first,a.second,b.first,b.second});
        }
    }
    std::set<int> requested(requests.begin(),requests.end());
    for(int i:requested) if(i<0 || i>=n) throw std::runtime_error("Refinement request out of range");
    std::vector<std::vector<unsigned char>> domain(n,std::vector<unsigned char>(no,1));
    std::vector<std::vector<std::array<int,6>>> signatures(n,std::vector<std::array<int,6>>(no));
    std::map<Signature,int> sig_ids;
    // Certify the parent first. An inverted cell is never repaired by merely
    // subdividing its inverted map; caller must supply a last-valid host.
    for(int i=0;i<n;++i) {
        const auto h=geometry(mesh,i);
        pillowgeom::Hex p;for(int k=0;k<8;++k)p[k]={h[k].x,h[k].y,h[k].z};
        if(!pillowgeom::valid(p,1)) throw std::runtime_error("Refinement requires valid hosts; parent "+std::to_string(i)+" is invalid or uncertifiable");
        for(int j=0;j<no;++j) {
            if(requested.count(i) && opt[j].pad!=144) domain[i][j]=0;
            for(int f=0;f<6;++f) {
                auto key=canonical_face(opt[j].faces[f],mesh.cells[i],f);
                auto r=sig_ids.emplace(key,int(sig_ids.size()));
                signatures[i][j][f]=r.first->second;
            }
        }
    }
    // Arc consistency plus bounded backtracking. Exact quadrilateral tilings
    // are the constraints, so no edge-only closure can leave hanging nodes.
    std::vector<std::pair<int,int>> trail;
    auto erase=[&](int i,int j){if(domain[i][j]){domain[i][j]=0;trail.push_back({i,j});}};
    std::size_t steps=0;
    const auto start=std::chrono::steady_clock::now();
    auto budget=[&](){
        if(++steps>max_steps || std::chrono::steady_clock::now()-start>std::chrono::seconds(7200))
            throw std::runtime_error("Refinement conformity search reached its work/time limit; mesh unchanged");
    };
    auto consistent=[&](){
        bool changed=true;
        while(changed) {
            budget();changed=false;
            for(const auto &l:links) for(int direction=0;direction<2;++direction) {
                int a=direction?l.b:l.a,fa=direction?l.fb:l.fa;
                int b=direction?l.a:l.b,fb=direction?l.fa:l.fb;
                std::set<int> supported;
                for(int j=0;j<no;++j)if(domain[b][j])supported.insert(signatures[b][j][fb]);
                for(int j=0;j<no;++j)if(domain[a][j] && !supported.count(signatures[a][j][fa])){erase(a,j);changed=true;}
            }
            std::size_t lower=0;
            for(int i=0;i<n;++i) {
                int first=-1;for(int j=0;j<no;++j)if(domain[i][j]){first=j;break;}
                if(first<0)return false;
                lower+=opt[first].children;
            }
            if(lower>max_elements)return false;
        }
        return true;
    };
    std::vector<int> choice(n,-1);
    // Geometry-valid choices are checked lazily, rather than evaluating all
    // 135 options on every host. Rejected options participate in backtracking.
    std::vector<std::vector<int8_t>> geom_ok(n,std::vector<int8_t>(no,-1));
    std::function<bool()> solve=[&](){
        if(!consistent())return false;
        // Commit every single-option host at this level instead of one recursion
        // level each: depth per host overflowed the stack on ~20k-host meshes.
        std::vector<int> forced;
        auto undo=[&](){for(int i:forced)choice[i]=-1;return false;};
        int selected=-1,best=no+1;
        for(int i=0;i<n;++i) {
            if(choice[i]>=0)continue;
            int count=std::accumulate(domain[i].begin(),domain[i].end(),0);
            if(count==1) {
                const int j=int(std::find(domain[i].begin(),domain[i].end(),1)-domain[i].begin());
                if(geom_ok[i][j]<0) geom_ok[i][j]=opt[j].pad==0 || validate_template_pad(opt[j].pad,geometry(mesh,i)).valid;
                if(!geom_ok[i][j])return undo();
                choice[i]=j;forced.push_back(i);
            } else if(count<best){selected=i;best=count;}
        }
        if(selected<0)return true;
        const int i=selected;
        for(int j=0;j<no;++j) if(domain[i][j]) {
            budget();
            if(geom_ok[i][j]<0) geom_ok[i][j]=opt[j].pad==0 || validate_template_pad(opt[j].pad,geometry(mesh,i)).valid;
            if(!geom_ok[i][j])continue;
            const auto checkpoint=trail.size();
            choice[i]=j;
            for(int k=0;k<no;++k)if(k!=j)erase(i,k);
            if(solve())return true;
            choice[i]=-1;
            while(trail.size()>checkpoint){auto p=trail.back();trail.pop_back();domain[p.first][p.second]=1;}
        }
        return undo();
    };
    if(!solve())throw std::runtime_error("No valid conforming template assignment fits the element budget; mesh unchanged");
    RefinedMesh result;
    result.mesh.points=mesh.points;result.requested=requested.size();result.pads.resize(n);
    std::map<Weights,int> nodes;
    for(int i=0;i<int(mesh.points.size());++i) {
        Weights w={{i,DEN}};nodes[w]=i;result.node_weights.push_back(w);
    }
    for(int i=0;i<n;++i) {
        const auto &o=opt[choice[i]];result.pads[i]=o.pad;result.replaced+=o.pad!=0;
        const auto h=geometry(mesh,i);
        for(const auto &child:o.cells) {
            HexConnectivity c;
            for(int k=0;k<8;++k) {
                auto key=weights(mesh.cells[i],child[k]);
                auto r=nodes.emplace(key,int(result.mesh.points.size()));
                if(r.second) {
                    result.mesh.points.push_back(trilinear_map(child[k],h));
                    result.node_weights.push_back(key);
                }
                c[k]=r.first->second;
            }
            result.mesh.cells.push_back(c);result.parent.push_back(i);
        }
    }
    // Independent final face incidence gate: single-owner faces must be on an
    // original exterior face; reject cracks and incompatible internal rings.
    std::map<Quad,std::vector<Quad>> final_faces;
    for(const auto &c:result.mesh.cells)for(const auto &f:FACES){
        Quad q;for(int k=0;k<4;++k)q[k]=c[f[k]];
        auto key=q;std::sort(key.begin(),key.end());final_faces[key].push_back(q);
    }
    for(const auto &p:final_faces) {
        if(p.second.size()>2)throw std::runtime_error("Refinement generated a non-manifold face");
        if(p.second.size()==2) {
            Quad reverse=p.second[1];std::reverse(reverse.begin(),reverse.end());
            auto cyclic=[](Quad q){auto m=std::min_element(q.begin(),q.end())-q.begin();Quad r;for(int k=0;k<4;++k)r[k]=q[(m+k)%4];return r;};
            if(cyclic(p.second[0])!=cyclic(reverse))throw std::runtime_error("Refinement face orientations disagree");
        } else {
            std::set<int> support;
            for(int id:p.first)for(auto w:result.node_weights[id])support.insert(w.first);
            bool exterior=false;
            // A full boundary quad has all four parent-face corners in its
            // combined support; reject unmatched faces inside the host mesh.
            if(support.size()==4) {
                Quad q;std::copy(support.begin(),support.end(),q.begin());
                auto it=owners.find(q);exterior=it!=owners.end() && it->second.size()==1;
            }
            if(!exterior)throw std::runtime_error("Refinement left an unmatched interior face");
        }
    }
    return result;
}
} // namespace refinement
