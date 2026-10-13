#include "hexa.h"
#include "pillow_geometry.h"
#include <mpi.h>
#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct Face {
    std::array<int,4> nodes;
    int element;
};
octant_t &element(hexa_tree_t *mesh, int id) {
    return *static_cast<octant_t *>(sc_array_index(&mesh->elements,id));
}
pillowgeom::Hex geometry(const octant_t &e, const std::vector<double> &coords) {
    pillowgeom::Hex p;
    for (int k=0; k<8; ++k) {
        int id=e.nodes[mgeom::H5_ORD[k]].id;
        for (int d=0; d<3; ++d) p[k][d]=coords[3*id+d];
    }
    return p;
}

// A failure on any rank cancels the insertion on every rank before mutation.
void check_all_ranks(hexa_tree_t *mesh, const std::string &error) {
    int failed=!error.empty(), any=failed;
    if (mesh->mpi_size>1) MPI_Allreduce(&failed,&any,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    if (any) throw std::runtime_error(error.empty()
        ? "Pillow cancelled: another MPI rank could not validate its insertion."
        : "Pillow cancelled before changing the mesh: "+error);
}

// Exterior faces of a rank can be material interfaces too. Match them using
// global lattice corners, since node IDs are local to each MPI rank. The pillow
// leaves these faces and their shared vertices intact; no buffer synchronization
// is needed after insertion.
void select_partition_interfaces(hexa_tree_t *mesh, const std::vector<Face> &boundary,
                                 std::vector<unsigned char> &selected,
                                 std::vector<unsigned char> &interface_node) {
    if (mesh->mpi_size<2) return;
    constexpr int width=16; // lattice coordinates, material, rank, element, boundary face
    std::vector<int> send;
    for (size_t fi=0; fi<boundary.size(); ++fi) {
        const auto &face=boundary[fi];
        std::array<std::array<int,3>,4> corners;
        for (int k=0; k<4; ++k) {
            auto *n=static_cast<octant_node_t *>(sc_array_index(&mesh->nodes,face.nodes[k]));
            corners[k]={n->x,n->y,n->z};
        }
        std::sort(corners.begin(),corners.end());
        for (const auto &c : corners) send.insert(send.end(),c.begin(),c.end());
        send.push_back(element(mesh,face.element).n_mat);
        send.push_back(mesh->mpi_rank);
        send.push_back(face.element);
        send.push_back(static_cast<int>(fi));
    }
    int count=static_cast<int>(send.size());
    std::vector<int> counts(mesh->mpi_size), offsets(mesh->mpi_size);
    MPI_Allgather(&count,1,MPI_INT,counts.data(),1,MPI_INT,MPI_COMM_WORLD);
    int total=0;
    for (int r=0; r<mesh->mpi_size; ++r) { offsets[r]=total; total+=counts[r]; }
    std::vector<int> received(total), order(total/width);
    MPI_Allgatherv(send.data(),count,MPI_INT,received.data(),counts.data(),offsets.data(),
                   MPI_INT,MPI_COMM_WORLD);
    for (int i=0; i<static_cast<int>(order.size()); ++i) order[i]=width*i;
    auto less=[&](int a,int b) {
        return std::lexicographical_compare(received.begin()+a,received.begin()+a+12,
                                             received.begin()+b,received.begin()+b+12);
    };
    std::sort(order.begin(),order.end(),less);
    for (size_t i=0; i<order.size();) {
        size_t j=i+1;
        while (j<order.size() && !less(order[i],order[j])) ++j;
        for (size_t a=i; a<j; ++a) if (received[order[a]+13]==mesh->mpi_rank)
            for (size_t b=i; b<j; ++b)
                if (received[order[a]+12]!=received[order[b]+12]) {
                    selected[received[order[a]+14]]=1;
                    for (int id : boundary[received[order[a]+15]].nodes) interface_node[id]=1;
                }
        i=j;
    }
}
} // namespace

void ApplyDoublePillowing(hexa_tree_t *mesh, std::vector<double> &coords,
                         std::vector<int> &nodes_b_mat) {
    if (!mesh) return;
    const size_t old_nodes=mesh->nodes.elem_count, old_elements=mesh->elements.elem_count;
    std::string error;
    if (coords.size()!=3*old_nodes) error="coordinate and node counts disagree";
    if (old_nodes>static_cast<size_t>(std::numeric_limits<int>::max()) ||
        old_elements>static_cast<size_t>(std::numeric_limits<int>::max()/6))
        error="mesh exceeds the supported 32-bit index range";
    check_all_ranks(mesh,error);
    printf("    Applying closed-cell double pillowing (validated before insertion)...\n");

    std::vector<Face> faces;
    faces.reserve(6*old_elements);
    double volume=0;
    for (size_t i=0; i<old_elements; ++i) {
        const octant_t &e=element(mesh,static_cast<int>(i));
        std::array<int,8> ids;
        for (int k=0; k<8; ++k) {
            ids[k]=e.nodes[k].id;
            if (ids[k]<0 || static_cast<size_t>(ids[k])>=old_nodes)
                error="out-of-range node in element "+std::to_string(i);
        }
        std::sort(ids.begin(),ids.end());
        if (std::adjacent_find(ids.begin(),ids.end())!=ids.end())
            error="repeated corner in element "+std::to_string(i);
        if (!error.empty()) break;
        auto p=geometry(e,coords);
        volume+=pillowgeom::jacobian(p,0.5,0.5,0.5);
        for (const auto &f : pillowgeom::faces) {
            Face face;
            face.element=static_cast<int>(i);
            for (int k=0; k<4; ++k) face.nodes[k]=e.nodes[mgeom::H5_ORD[f[k]]].id;
            std::sort(face.nodes.begin(),face.nodes.end());
            faces.push_back(face);
        }
    }
    check_all_ranks(mesh,error);
    if (mesh->mpi_size>1) {
        double local=volume;
        MPI_Allreduce(&local,&volume,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    }
    const int orientation=volume>=0 ? 1 : -1;
    std::sort(faces.begin(),faces.end(),[](const Face &a,const Face &b) { return a.nodes<b.nodes; });
    std::vector<unsigned char> selected(old_elements,0), interface_node(old_nodes,0);
    std::vector<Face> boundary;
    std::set<std::array<int,4>> interface_faces;   // sorted node ids of material-interface faces
    size_t interfaces=0;
    for (size_t i=0; i<faces.size();) {
        size_t j=i+1;
        while (j<faces.size() && faces[j].nodes==faces[i].nodes) ++j;
        if (j-i>2) { error="non-manifold face at element "+std::to_string(faces[i].element); break; }
        if (j-i==1) boundary.push_back(faces[i]);
        if (j-i==2 && element(mesh,faces[i].element).n_mat!=element(mesh,faces[i+1].element).n_mat) {
            selected[faces[i].element]=selected[faces[i+1].element]=1;
            interface_faces.insert(faces[i].nodes);
            for (int id : faces[i].nodes) interface_node[id]=1;
            ++interfaces;
        }
        i=j;
    }
    check_all_ranks(mesh,error);
    std::vector<Face>().swap(faces);
    select_partition_interfaces(mesh,boundary,selected,interface_node);
    std::vector<Face>().swap(boundary);

    struct Insertion { int host; pillowgeom::Hex inner; std::array<int,6> rotations; double t[6]; };
    double min_inset=0.175;
    std::vector<Insertion> insertions;
    insertions.reserve(std::count(selected.begin(),selected.end(),1));
    for (size_t i=0; i<old_elements; ++i) {
        auto p=geometry(element(mesh,static_cast<int>(i)),coords);
        if (!pillowgeom::valid(p,orientation)) {
            error="input element "+std::to_string(i)+" has a non-positive or uncertifiable Jacobian";
            break;
        }
        if (!selected[i]) continue;
        Insertion insertion;
        insertion.host=static_cast<int>(i);
        double inset;
        const octant_t &host=element(mesh,static_cast<int>(i));
        // pillowFaceInset = 1 (default): half-cell shell on interface faces (default_insets);
        // 0: the old uniform 0.175 box
        const bool half=mesh->input.pillowFaceInset || mesh->input.pillowFirst;
        unsigned mask=0x3F;
        if (half) mask=0;
        for (int f=0; f<6 && mesh->input.pillowFaceInset; ++f) {
            std::array<int,4> key;
            for (int k=0; k<4; ++k) key[k]=host.nodes[mgeom::H5_ORD[pillowgeom::faces[f][k]]].id;
            std::sort(key.begin(),key.end());
            if (interface_faces.count(key)) mask|=1u<<f;
        }
        static const bool fixed_inset=getenv("PILLOW_FIXED_INSET")!=nullptr;   // A/B: old uniform backoff
        pillowgeom::default_insets(mask,half,insertion.t);   // also the lattice layout below
        if (!(fixed_inset ? pillowgeom::build(p,orientation,insertion.inner,&inset,mask,half)
                          : pillowgeom::build_optimized(p,orientation,insertion.inner,&inset,mask,half))) {
            error="no valid closed pillow fits inside element "+std::to_string(i);
            break;
        }
        min_inset=std::min(min_inset,inset);
        for (int f=0; f<6; ++f)
            insertion.rotations[f]=pillowgeom::shell_rotation(p,insertion.inner,f,orientation);
        insertions.push_back(insertion);
    }
    check_all_ranks(mesh,error);
    if (getenv("PILLOW_DT_DIAG")) {
        // which hosts give the thinnest children, and is the host itself already thin there
        struct Row { double child, host; int host_id, which; };
        std::vector<Row> rows;
        for (const auto &ins : insertions) {
            auto p=geometry(element(mesh,ins.host),coords);
            Row r{pillowgeom::size_metric(ins.inner),pillowgeom::size_metric(p),ins.host,6};
            for (int f=0; f<6; ++f) {
                double s=pillowgeom::size_metric(pillowgeom::shell_geometry(p,ins.inner,f,ins.rotations[f]));
                if (s<r.child) { r.child=s; r.which=f; }
            }
            rows.push_back(r);
        }
        std::sort(rows.begin(),rows.end(),[](const Row &a,const Row &b) { return a.child<b.child; });
        printf("    PILLOW_DT_DIAG: worst hosts by child size metric (sigma_min of J, m)\n");
        for (size_t i=0; i<std::min<size_t>(12,rows.size()); ++i)
            printf("      host %7d  child %-5s %.4e   host %.4e   ratio %.3f\n",rows[i].host_id,
                   rows[i].which==6 ? "core" : ("f"+std::to_string(rows[i].which)).c_str(),
                   rows[i].child,rows[i].host,rows[i].child/rows[i].host);
    }
    // Every child has a positive Jacobian throughout its reference cube. The
    // seven children share their internal faces and preserve the host boundary.
    // Their positive maps make this a conforming subdivision of the host.
    long long additions[2]={static_cast<long long>(8*insertions.size()),
                            static_cast<long long>(6*insertions.size())};
    long long global_additions[2]={additions[0],additions[1]};
    if (mesh->mpi_size>1)
        MPI_Allreduce(additions,global_additions,2,MPI_LONG_LONG_INT,MPI_SUM,MPI_COMM_WORLD);
    if (old_nodes+additions[0]>static_cast<size_t>(std::numeric_limits<int>::max()) ||
        old_elements+additions[1]>static_cast<size_t>(std::numeric_limits<int>::max()))
        error="pillow would exceed the supported 32-bit index range";
    check_all_ranks(mesh,error);
    if (global_additions[0]==0) {
        printf("    Pillow: no material interfaces; %zu elements validated, 0 invalid.\n",old_elements);
        return;
    }

    constexpr int factor=4;
    const int bound=std::numeric_limits<int>::max()/factor;
    auto fits=[&](int x) { return x>=-bound && x<=bound; };
    for (size_t id=0; id<old_nodes; ++id) {
        auto *n=static_cast<octant_node_t *>(sc_array_index(&mesh->nodes,id));
        if (!fits(n->x) || !fits(n->y) || !fits(n->z)) error="integer node coordinates would overflow";
    }
    for (size_t i=0; i<old_elements; ++i) {
        const auto &e=element(mesh,static_cast<int>(i));
        if (!fits(e.x) || !fits(e.y) || !fits(e.z)) error="integer element coordinates would overflow";
    }
    if (!fits(mesh->x_start) || !fits(mesh->x_end) || !fits(mesh->y_start) || !fits(mesh->y_end))
        error="integer partition bounds would overflow";
    check_all_ranks(mesh,error);

    // Reserve before altering connectivity. Keep the ownership of existing
    // shared nodes; newly created nodes are strictly inside a local host.
    coords.reserve(3*(old_nodes+additions[0]));
    std::vector<int> registered=nodes_b_mat;
    std::vector<unsigned char> marked=interface_node;
    for (int id : registered) if (id>=0 && static_cast<size_t>(id)<old_nodes) marked[id]=0;
    for (size_t id=0; id<old_nodes; ++id) if (marked[id]) registered.push_back(static_cast<int>(id));
    auto *owners=static_cast<int32_t *>(std::malloc((old_nodes+additions[0])*sizeof(int32_t)));
    if (!owners && old_nodes+additions[0]) error="cannot allocate pillow node ownership";
    try { check_all_ranks(mesh,error); }
    catch (...) { std::free(owners); throw; }
    for (size_t id=0; id<old_nodes+additions[0]; ++id)
        owners[id]=(id<old_nodes && mesh->part_nodes) ? mesh->part_nodes[id] : mesh->mpi_rank;

    // Integer coordinates are topology metadata, not physical positions. Four
    // slots per old lattice unit leave strictly interior slots for the new nodes.
    for (size_t id=0; id<old_nodes; ++id) {
        auto *n=static_cast<octant_node_t *>(sc_array_index(&mesh->nodes,id));
        n->x*=factor; n->y*=factor; n->z*=factor;
    }
    for (size_t i=0; i<old_elements; ++i) {
        auto &e=element(mesh,static_cast<int>(i));
        e.x*=factor; e.y*=factor; e.z*=factor;
        for (auto &n : e.nodes) { n.x*=factor; n.y*=factor; n.z*=factor; }
    }
    mesh->x_start*=factor; mesh->x_end*=factor;
    mesh->y_start*=factor; mesh->y_end*=factor;
    for (const auto &insertion : insertions) {
        // Copy before sc_array_push: a growth of either array may relocate it.
        const octant_t host=element(mesh,insertion.host);
        std::array<octant_node_t,8> inner_nodes;
        // lattice position of the core: the host's lattice box restricted by the default insets
        // (not the optimized ones), so the reference mesh shows the layout -- e.g. the half-cell
        // shell on an interface face spans 2 of the host's 4 slots
        pillowgeom::Hex lat, lat_core;
        for (int k=0; k<8; ++k) {
            const auto &o=host.nodes[mgeom::H5_ORD[k]];
            lat[k]={static_cast<double>(o.x),static_cast<double>(o.y),static_cast<double>(o.z)};
        }
        pillowgeom::core_from_insets(lat,insertion.t,lat_core);
        for (int k=0; k<8; ++k) {
            octant_node_t node{};
            node.id=static_cast<int>(mesh->nodes.elem_count);
            node.x=static_cast<int>(std::lround(lat_core[k][0]));
            node.y=static_cast<int>(std::lround(lat_core[k][1]));
            node.z=static_cast<int>(std::lround(lat_core[k][2]));
            *static_cast<octant_node_t *>(sc_array_push(&mesh->nodes))=node;
            inner_nodes[k]=node;
            coords.insert(coords.end(),insertion.inner[k].begin(),insertion.inner[k].end());
        }
        pillow_host_t rec;
        for (int k=0; k<8; ++k) { rec.outer[k]=host.nodes[mgeom::H5_ORD[k]].id; rec.core[k]=inner_nodes[k].id; }
        std::copy(insertion.t,insertion.t+6,rec.t);
        mesh->pillow_hosts.push_back(rec);
        auto &core=element(mesh,insertion.host);
        for (int k=0; k<8; ++k) core.nodes[mgeom::H5_ORD[k]]=inner_nodes[k];
        for (auto &s : core.surf) s.ext=false;
        core.boundary=false;
        for (int f=0; f<6; ++f) {
            octant_t shell{};
            shell.id=static_cast<int64_t>(mesh->elements.elem_count);
            shell.n_mat=host.n_mat; shell.level=host.level; shell.father=host.father;
            shell.x=host.x; shell.y=host.y; shell.z=host.z;
            shell.pad=host.pad; shell.pml_id=host.pml_id;
            const int shift=insertion.rotations[f];
            const auto corners=pillowgeom::shell_corners(f,shift);
            for (int k=0; k<8; ++k) {
                const int j=corners[k];
                shell.nodes[mgeom::H5_ORD[k]]=j<8 ? inner_nodes[j] : host.nodes[mgeom::H5_ORD[j-8]];
            }
            *static_cast<octant_t *>(sc_array_push(&mesh->elements))=shell;
        }
    }
    std::free(mesh->part_nodes);
    mesh->part_nodes=owners;
    nodes_b_mat.swap(registered);
    mesh->local_n_nodes=static_cast<int32_t>(mesh->nodes.elem_count);
    mesh->local_n_elements=static_cast<int32_t>(mesh->elements.elem_count);
    mesh->total_n_nodes+=global_additions[0];
    mesh->total_n_elements+=global_additions[1];
    printf("    Pillow minimum reference inset: %.6g\n",min_inset);
    printf("    Pillow interfaces: %zu local faces, %zu host cells\n",interfaces,insertions.size());
    printf("    Pillowing summary: %lld new nodes created, %lld new pillow elements added\n",additions[0],additions[1]);
    printf("    Pillow validation: %zu elements, 0 invalid (positive Jacobian throughout each cell)\n",
           mesh->elements.elem_count);
}
