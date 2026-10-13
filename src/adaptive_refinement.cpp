#include "adaptive_refinement.h"
#include "conforming_refinement.h"
#include "pillow_geometry.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>

namespace {
constexpr int DEN=24*24*24;
double min_sj(const refinement::HexCorners &h) {
    double x[8],y[8],z[8];
    for(int k=0;k<8;++k){x[k]=h[k].x;y[k]=h[k].y;z[k]=h[k].z;}
    return mgeom::hex_min_corner_sj(x,y,z);
}
refinement::HexCorners get_hex(const refinement::MeshData &m,int i) {
    refinement::HexCorners h;for(int k=0;k<8;++k)h[k]=m.points[m.cells[i][k]];return h;
}

// Refine-then-reproject: children were interpolated from the projected parent,
// so new nodes on a material interface sit on the parent's bilinear face, not
// on the surface. Snap them vertically onto the nearest interface sheet, spread
// the dz into the other new nodes (harmonic extension, old nodes fixed), and
// halve dz on the nodes of any child that loses full-Jacobian validity.
// Only z moves, and only new nodes (id >= nn); the sea-lid rules match
// ProjectFreeNodes (z >= 0 stays, otherwise z <= -1 m).
void reproject_new_interface_nodes(hexa_tree_t *mesh,const std::vector<std::array<int,8>> &cells,int nn,
                                   const std::vector<int> &iface,std::vector<double> &X) {
    const int np=int(X.size()/3);
    std::vector<double> dz(np,0.0);
    std::vector<char> fixed(np,1);
    for(int i=nn;i<np;++i)fixed[i]=0;
    std::vector<int> surfaces;
    for(int k=0;k<int(std::max<std::size_t>(1,mesh->gdata_vec.size()));++k)surfaces.push_back(k);
    int snapped=0,missed=0;double maxdz=0,sumdz=0;
    for(int i:iface) {
        if(i<nn)continue;
        fixed[i]=1;
        const double x=X[3*i],y=X[3*i+1],z0=X[3*i+2];
        if(z0>=-1e-6)continue;  // sea-surface lid node
        double best=0;bool found=false;
        for(int id:surfaces){double z;if(eval_gts_height_near(mesh,id,x,y,z0,z) && (!found || std::fabs(z-z0)<std::fabs(best-z0))){best=z;found=true;}}
        if(!found){++missed;continue;}
        best=std::min(best,-1.0);
        dz[i]=best-z0;++snapped;maxdz=std::max(maxdz,std::fabs(dz[i]));sumdz+=std::fabs(dz[i]);
    }
    std::vector<std::vector<int>> adj(np);
    static const int edges[12][2]={{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}};
    for(const auto &c:cells)for(const auto &e:edges){adj[c[e[0]]].push_back(c[e[1]]);adj[c[e[1]]].push_back(c[e[0]]);}
    // ponytail: fixed 50 Gauss-Seidel sweeps; correction only needs to reach across one refined host (3-4 nodes)
    for(int it=0;it<50;++it)for(int i=nn;i<np;++i)if(!fixed[i] && !adj[i].empty()){
        double a=0;for(int j:adj[i])a+=dz[j];dz[i]=a/adj[i].size();
    }
    auto hex=[&](const std::array<int,8> &c){pillowgeom::Hex h;for(int k=0;k<8;++k)h[k]={X[3*c[k]],X[3*c[k]+1],X[3*c[k]+2]+dz[c[k]]};return h;};
    int retreats=0,rounds=0,bad=0;
    for(;rounds<12;++rounds) {
        std::vector<char> shrink(np,0);bad=0;
        for(const auto &c:cells) {
            bool moved=false;for(int k:c)moved|=dz[k]!=0;
            if(!moved || pillowgeom::valid(hex(c),1))continue;
            ++bad;for(int k:c)shrink[k]=1;
        }
        if(!bad)break;
        for(int i=nn;i<np;++i)if(shrink[i] && dz[i]!=0){dz[i]=rounds<11?0.5*dz[i]:0.0;++retreats;}
    }
    for(int i=nn;i<np;++i)X[3*i+2]+=dz[i];
    printf("    Refinement reprojection: %d new interface nodes snapped (%d missed), mean |dz| %.1f m, max |dz| %.1f m; %d node retreats over %d rounds, %d cells still invalid\n",
           snapped,missed,snapped?sumdz/snapped:0.0,maxdz,retreats,rounds,bad);
}
}

void ApplyAdaptiveRefinement(hexa_tree_t *mesh,std::vector<double> &coords,
                             std::vector<int> &nodes_b_mat) {
    if(!mesh->input.refinement)return;
    if(mesh->mpi_size!=1)throw std::runtime_error("refinement=1 currently requires one MPI rank");
    if(mesh->input.PML)throw std::runtime_error("refinement=1 currently requires PML=0; refined PML is not qualified");
    refinement::MeshData data;
    const int nn=mesh->nodes.elem_count,ne=mesh->elements.elem_count;
    if(coords.size()!=3*std::size_t(nn))throw std::runtime_error("Refinement coordinate count mismatch");
    for(int i=0;i<nn;++i)data.points.push_back({coords[3*i],coords[3*i+1],coords[3*i+2]});
    std::vector<octant_t> parents;
    std::set<int> request(mesh->refinement_projection_requests.begin(),mesh->refinement_projection_requests.end());
    double before=1;
    for(int i=0;i<ne;++i) {
        auto e=*static_cast<octant_t *>(sc_array_index(&mesh->elements,i));parents.push_back(e);
        refinement::HexConnectivity c;for(int k=0;k<8;++k)c[k]=e.nodes[mgeom::H5_ORD[k]].id;
        data.cells.push_back(c);
        const double q=min_sj(get_hex(data,i));before=std::min(before,q);
        if(!std::isfinite(q) || q<mesh->input.refinementMinSJ)request.insert(i);
    }
    printf("    Adaptive refinement: %zu requested hosts (raw projection failures or corner SJ < %.6g), %d hosts, minimum SJ %.6g\n",
           request.size(),mesh->input.refinementMinSJ,ne,before);
    if(request.empty()) {printf("    Adaptive refinement: no poor hosts, mesh unchanged\n");return;}
    auto result=refinement::refine_conforming(data,{request.begin(),request.end()},mesh->input.refinementMaxElements,100000000);
    // Account for the actual subsequent closed-pillow insertion before mutation.
    std::map<std::array<int,4>,int> faces;
    std::set<int> pillow_hosts;
    for(int i=0;i<int(result.mesh.cells.size());++i) for(const auto &f:pillowgeom::faces) {
        std::array<int,4> key;for(int k=0;k<4;++k)key[k]=result.mesh.cells[i][f[k]];
        std::sort(key.begin(),key.end());auto r=faces.emplace(key,i);
        if(!r.second && parents[result.parent[i]].n_mat!=parents[result.parent[r.first->second]].n_mat) {
            pillow_hosts.insert(i);pillow_hosts.insert(r.first->second);
        }
    }
    const auto final_size=result.mesh.cells.size()+(mesh->input.movingNodes?6*pillow_hosts.size():0);
    if(final_size>std::size_t(mesh->input.refinementMaxElements))
        throw std::runtime_error("Refinement plus closed pillow exceeds refinementMaxElements; mesh unchanged");
    // Exact integer logical interpolation. This scale makes every weighted
    // logical coordinate integral and preserves boundary planes. It is checked
    // for overflow before any live storage changes.
    std::vector<octant_node_t> nodes(result.mesh.points.size());
    auto scaled=[](int value){
        int64_t out=int64_t(value)*DEN;
        if(out<std::numeric_limits<int32_t>::min() || out>std::numeric_limits<int32_t>::max())
            throw std::runtime_error("Refinement logical coordinate overflow; mesh unchanged");
        return int(out);
    };
    for(std::size_t i=0;i<nodes.size();++i) {
        auto &node=nodes[i];node={};node.id=i;node.color=-1;node.fixed=0;
        int64_t xyz[3]={0,0,0};
        for(auto w:result.node_weights[i]) {
            const auto &old=*static_cast<octant_node_t *>(sc_array_index(&mesh->nodes,w.first));
            xyz[0]+=int64_t(w.second)*old.x;xyz[1]+=int64_t(w.second)*old.y;xyz[2]+=int64_t(w.second)*old.z;
        }
        for(int d=0;d<3;++d)if(xyz[d]<std::numeric_limits<int32_t>::min() || xyz[d]>std::numeric_limits<int32_t>::max())
            throw std::runtime_error("Refinement logical coordinate overflow; mesh unchanged");
        node.x=xyz[0];node.y=xyz[1];node.z=xyz[2];
        if(i<std::size_t(nn)) {const auto old=*static_cast<octant_node_t *>(sc_array_index(&mesh->nodes,i));node.color=old.color;node.fixed=old.fixed;}
    }
    const int xs=scaled(mesh->x_start),xe=scaled(mesh->x_end),ys=scaled(mesh->y_start),ye=scaled(mesh->y_end);
    const int nx=scaled(mesh->ncellx),ny=scaled(mesh->ncelly),nz=scaled(mesh->ncellz),mz=scaled(mesh->max_z),step=scaled(mesh->max_step);
    std::vector<octant_t> elements;
    double after=1;
    for(std::size_t i=0;i<result.mesh.cells.size();++i) {
        const int pi=result.parent[i];const auto &parent=parents[pi];
        octant_t e=parent;e.id=i;e.father=pi;
        e.level=parent.level+(result.pads[pi]!=0);e.tem=result.pads[pi]?refinement::get_pad_orientation(result.pads[pi])->template_id:0;
        e.pad=0;e.x=scaled(parent.x);e.y=scaled(parent.y);e.z=scaled(parent.z);
        for(auto &edge:e.edge) edge={};
        for(auto &face:e.surf)face={};e.boundary=false;
        for(int k=0;k<8;++k)e.nodes[mgeom::H5_ORD[k]]=nodes[result.mesh.cells[i][k]];
        // Preserve physical boundary tags on every descendant of an old face.
        for(int f=0;f<6;++f) {
            std::set<int> support;
            for(int k=0;k<4;++k)for(auto w:result.node_weights[e.nodes[FaceNodesMap[f][k]].id])support.insert(w.first);
            for(int pf=0;pf<6;++pf) {
                std::set<int> oldface;for(int k=0;k<4;++k)oldface.insert(parent.nodes[FaceNodesMap[pf][k]].id);
                if(support==oldface){e.surf[f]=parent.surf[pf];e.boundary|=e.surf[f].ext;}
            }
        }
        elements.push_back(e);after=std::min(after,min_sj(get_hex(result.mesh,i)));
    }
    // Prepare all allocations and coordinates before commit.
    std::vector<double> candidate;candidate.reserve(3*nodes.size());
    for(auto p:result.mesh.points)candidate.insert(candidate.end(),{p.x,p.y,p.z});
    std::vector<int> interface_nodes;
    std::map<std::array<int,4>,int> iface;
    for(int i=0;i<int(elements.size());++i)for(const auto &f:pillowgeom::faces) {
        std::array<int,4> key;for(int k=0;k<4;++k)key[k]=result.mesh.cells[i][f[k]];
        std::sort(key.begin(),key.end());auto r=iface.emplace(key,i);
        if(!r.second && elements[i].n_mat!=elements[r.first->second].n_mat)interface_nodes.insert(interface_nodes.end(),key.begin(),key.end());
    }
    std::sort(interface_nodes.begin(),interface_nodes.end());interface_nodes.erase(std::unique(interface_nodes.begin(),interface_nodes.end()),interface_nodes.end());
    {
        const auto &cells=result.mesh.cells;
        reproject_new_interface_nodes(mesh,cells,nn,interface_nodes,candidate);
        after=1;
        for(const auto &c:cells){double x[8],y[8],z[8];for(int k=0;k<8;++k){x[k]=candidate[3*c[k]];y[k]=candidate[3*c[k]+1];z[k]=candidate[3*c[k]+2];}after=std::min(after,mgeom::hex_min_corner_sj(x,y,z));}
    }
    auto *global=static_cast<int64_t *>(malloc(nodes.size()*sizeof(int64_t)));
    auto *partition=static_cast<int32_t *>(malloc(nodes.size()*sizeof(int32_t)));
    if(!global || !partition){free(global);free(partition);throw std::runtime_error("Refinement allocation failed");}
    for(std::size_t i=0;i<nodes.size();++i){global[i]=i;partition[i]=0;}
    sc_array_resize(&mesh->nodes,nodes.size());sc_array_resize(&mesh->elements,elements.size());
    for(std::size_t i=0;i<nodes.size();++i)*static_cast<octant_node_t *>(sc_array_index(&mesh->nodes,i))=nodes[i];
    for(std::size_t i=0;i<elements.size();++i)*static_cast<octant_t *>(sc_array_index(&mesh->elements,i))=elements[i];
    free(mesh->global_id);mesh->global_id=global;free(mesh->part_nodes);mesh->part_nodes=partition;
    // Old grouping and edge identities must not be used after replacement.
    if(mesh->oct.elem_size)sc_array_reset(&mesh->oct);
    coords.swap(candidate);nodes_b_mat.swap(interface_nodes);
    mesh->local_n_nodes=mesh->total_n_nodes=nodes.size();mesh->local_n_elements=mesh->total_n_elements=elements.size();
    mesh->x_start=xs;mesh->x_end=xe;mesh->y_start=ys;mesh->y_end=ye;
    mesh->ncellx=nx;mesh->ncelly=ny;mesh->ncellz=nz;mesh->max_z=mz;mesh->max_step=step;
    mesh->refinement_projection_requests.clear();
    printf("    Adaptive refinement: %zu requested, %zu replaced including transitions; hosts %d -> %zu, nodes %d -> %zu, predicted post-pillow %zu\n",
           result.requested,result.replaced,ne,elements.size(),nn,nodes.size(),final_size);
    printf("    Adaptive refinement comparison: minimum corner SJ %.9g -> %.9g (positive full-volume Jacobians and conforming faces checked)\n",before,after);
}
