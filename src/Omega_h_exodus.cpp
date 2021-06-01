#include "Omega_h_mpi.h"

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#endif
#include <exodusII.h>
#ifdef __clang__
#pragma clang diagnostic pop
#endif

#include <algorithm>
#include <iostream>
#include <set>
#include <sstream>

#include "Omega_h_align.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_build.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mark.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_dbg.hpp"

namespace Omega_h {

#define CALL(f)                                                                \
  do {                                                                         \
    auto f_err = (f);                                                          \
    if (f_err != 0) {                                                          \
      const char* errmsg;                                                      \
      const char* errfunc;                                                     \
      int errnum;                                                              \
      ex_get_err(&errmsg, &errfunc, &errnum);                                  \
      Omega_h_fail("Exodus call %s failed (%d): %s: %s\n", #f, errnum,         \
          errfunc, errmsg);                                                    \
    }                                                                          \
  } while (0)

namespace exodus {

static void get_elem_type_info(
    std::string const& type, int* p_dim, Omega_h_Family* p_family) {
  if (type == "tri3") {
    *p_dim = 2;
    *p_family = OMEGA_H_SIMPLEX;
  } else if (type == "TRI") {
    *p_dim = 2;
    *p_family = OMEGA_H_SIMPLEX;
  } else if (type == "TRI3") {
    *p_dim = 2;
    *p_family = OMEGA_H_SIMPLEX;
  } else if (type == "tetra4") {
    *p_dim = 3;
    *p_family = OMEGA_H_SIMPLEX;
  } else if (type == "TETRA") {
    *p_dim = 3;
    *p_family = OMEGA_H_SIMPLEX;
  } else if (type == "TET4") {
    *p_dim = 3;
    *p_family = OMEGA_H_SIMPLEX;
  } else {
    Omega_h_fail("Unsupported Exodus element type \"%s\"\n", type.c_str());
  }
}

// subtracts one and maps from Exodus
// side ordering to Omega_h
static OMEGA_H_INLINE int side_exo2osh(
    Omega_h_Family family, int dim, int side) {
  switch (family) {
    case OMEGA_H_SIMPLEX:
      switch (dim) {
        case 2:
          // seeing files from CUBIT with triangle sides in {3,4,5}...
          // no clue what thats about, just modulo and move on
          return (side) % 3;
        case 3:
          switch (side) {
            case 1:
              return 1;
            case 2:
              return 2;
            case 3:
              return 3;
            case 4:
              return 0;
          }
      }
      return -1;
    case OMEGA_H_HYPERCUBE:
      return -1;  // needs to be filled in!
  }
  return -1;
}

// from Omega_h
// side ordering to Exodus and adds one
static OMEGA_H_INLINE int side_osh2exo(int dim, int side) {
  switch (dim) {
    case 2:
      switch (side) {
        case 0:
          return 1;
        case 1:
          return 2;
        case 2:
          return 3;
      }
      return -1;
    case 3:
      switch (side) {
        case 1:
          return 1;
        case 2:
          return 2;
        case 3:
          return 3;
        case 0:
          return 4;
      }
      return -1;
  }
  return -1;
}

int open(filesystem::path const& path, bool verbose) {
  auto comp_ws = int(sizeof(Real));
  int io_ws = 0;
  float version;
  auto mode = EX_READ | EX_MAPS_INT64_API;
  auto exodus_file = ex_open(path.c_str(), mode, &comp_ws, &io_ws, &version);
  if (exodus_file < 0)
    Omega_h_fail("can't open Exodus file %s\n", path.c_str());
  if (verbose) {
    std::cout << "ex_open(" << path << ")\n";
    std::cout << "  comp_ws: " << comp_ws << '\n';
    std::cout << "  io_ws: " << io_ws << '\n';
    std::cout << "  version: " << version << '\n';
  }
  return exodus_file;
}

void close(int exodus_file) { CALL(ex_close(exodus_file)); }

int get_num_time_steps(int exodus_file) {
  return int(ex_inquire_int(exodus_file, EX_INQ_TIME));
}

static void setup_names(
    int nnames, std::vector<char>& storage, std::vector<char*>& ptrs) {
  constexpr auto max_name_length = MAX_STR_LENGTH + 1;
  storage = std::vector<char>(std::size_t(nnames * max_name_length), '\0');
  ptrs = std::vector<char*>(std::size_t(nnames), nullptr);
  for (int i = 0; i < nnames; ++i) {
    ptrs[std::size_t(i)] = storage.data() + max_name_length * i;
  }
}

void read_nodal_fields(int exodus_file, Mesh* mesh, int time_step,
    std::string const& prefix, std::string const& postfix, bool verbose) {
  int num_nodal_vars;
  CALL(ex_get_variable_param(exodus_file, EX_NODAL, &num_nodal_vars));
  if (verbose) std::cout << "P" << mesh->comm()->rank() << ": " << num_nodal_vars << " nodal variables\n";
  if (num_nodal_vars == 0) return;
  std::vector<char> names_memory;
  std::vector<char*> name_ptrs;
  setup_names(num_nodal_vars, names_memory, name_ptrs);
  CALL(ex_get_variable_names(
      exodus_file, EX_NODAL, num_nodal_vars, name_ptrs.data()));
  for (int i = 0; i < num_nodal_vars; ++i) {
    auto name = name_ptrs[std::size_t(i)];
    if (verbose)
      std::cout << "P" << mesh->comm()->rank() << ": Loading nodal variable \"" << name << "\" at time step "
                << time_step << '\n';
    auto name_osh = prefix + std::string(name) + postfix;
    HostWrite<double> host_write(mesh->nverts(), name_osh);
    CALL(ex_get_var(exodus_file, time_step + 1, EX_NODAL, i + 1, /*obj_id*/ 0,
        mesh->nverts(), host_write.data()));
    auto device_write = host_write.write();
    auto device_read = Reals(device_write);
    mesh->add_tag(VERT, name_osh, 1, device_read);
  }
}

void read_mesh(int file, Mesh* mesh, bool verbose, int classify_with) {
  begin_code("exodus::read_mesh");
  ex_init_params init_params;
  CALL(ex_get_init_ext(file, &init_params));
  if (verbose) {
    std::ostringstream oss;
    oss << "init params:\n";
    oss << " Exodus ID " << file << '\n';
    oss << " Title " << init_params.title << '\n';
    oss << " num_dim " << init_params.num_dim << '\n';
    oss << " num_nodes " << init_params.num_nodes << '\n';
    oss << " num_elem " << init_params.num_elem << '\n';
    oss << " num_elem_blk " << init_params.num_elem_blk << '\n';
    oss << " num_node_sets " << init_params.num_node_sets << '\n';
    oss << " num_side_sets " << init_params.num_side_sets << '\n';
    std::cout << oss.str() << std::endl;
  }
  std::vector<int> block_ids(std::size_t(init_params.num_elem_blk));
  CALL(ex_get_ids(file, EX_ELEM_BLOCK, block_ids.data()));
  std::vector<char> block_names_memory;
  std::vector<char*> block_names;
  setup_names(int(init_params.num_elem_blk), block_names_memory, block_names);
  CALL(ex_get_names(file, EX_ELEM_BLOCK, block_names.data()));
  HostWrite<LO> h_conn;
  Write<LO> elem_class_ids_w(LO(init_params.num_elem));
  LO elem_start = 0;
  int family_int = -1;
  int dim = -1;
  for (size_t i = 0; i < block_ids.size(); ++i) {
    char elem_type[MAX_STR_LENGTH + 1];
    int nentries;
    int nnodes_per_entry;
    int nedges_per_entry;
    int nfaces_per_entry;
    int nattr_per_entry;
    CALL(ex_get_block(file, EX_ELEM_BLOCK, block_ids[i], elem_type, &nentries,
        &nnodes_per_entry, &nedges_per_entry, &nfaces_per_entry,
        &nattr_per_entry));
    if (verbose) {
      std::cout << "P" << mesh->comm()->rank() << ": block " << block_ids[i] << " \"" << block_names[i] << "\""
                << " has " << nentries << " elements of type " << elem_type
                << std::endl;
    }
    /* some pretty weird blocks from the CDFEM people... */
    if (std::string("NULL") == elem_type && nentries == 0) continue;
    int dim_from_type;
    Omega_h_Family family_from_type;
    get_elem_type_info(elem_type, &dim_from_type, &family_from_type);
    if (family_int == -1) family_int = family_from_type;
    OMEGA_H_CHECK_OP(family_int, ==, family_from_type);
    if (dim == -1) dim = dim_from_type;
    OMEGA_H_CHECK_OP(dim, ==, dim_from_type);
    auto deg = element_degree(Omega_h_Family(family_int), dim, VERT);
    OMEGA_H_CHECK_OP(nnodes_per_entry, ==, deg);
    if (!h_conn.exists())
      h_conn =
          decltype(h_conn)(LO(init_params.num_elem * deg), "host connectivity");
    if (nedges_per_entry < 0) nedges_per_entry = 0;
    if (nfaces_per_entry < 0) nfaces_per_entry = 0;
    std::vector<int> edge_conn(std::size_t(nentries * nedges_per_entry));
    std::vector<int> face_conn(std::size_t(nentries * nfaces_per_entry));
    CALL(ex_get_conn(file, EX_ELEM_BLOCK, block_ids[i],
        h_conn.data() + elem_start * nnodes_per_entry, edge_conn.data(),
        face_conn.data()));
    auto region_id = block_ids[i];
    auto f0 = OMEGA_H_LAMBDA(LO entry) {
      elem_class_ids_w[elem_start + entry] = region_id;
    };
    parallel_for(nentries, f0, "set_elem_class_ids");
    mesh->class_sets[block_names[i]].push_back({I8(dim), region_id});
    elem_start += nentries;
  }
  OMEGA_H_CHECK_OP(elem_start, ==, init_params.num_elem);
  Omega_h_Family family = Omega_h_Family(family_int);
  auto conn = subtract_from_each(LOs(h_conn.write()), 1);
  HostWrite<Real> h_coord_blk[3];
  for (Int i = 0; i < dim; ++i) {
    h_coord_blk[i] = HostWrite<Real>(LO(init_params.num_nodes));
  }
  CALL(ex_get_coord(file, h_coord_blk[0].data(), h_coord_blk[1].data(),
      h_coord_blk[2].data()));
  HostWrite<Real> h_coords(LO(init_params.num_nodes * dim));
  for (LO i = 0; i < init_params.num_nodes; ++i) {
    for (Int j = 0; j < dim; ++j) {
      h_coords[i * dim + j] = h_coord_blk[j][i];
    }
  }
  auto coords = Reals(h_coords.write());
  build_from_elems_and_coords(mesh, OMEGA_H_SIMPLEX, dim, conn, coords);
  classify_elements(mesh);
  std::vector<int> side_set_ids(std::size_t(init_params.num_side_sets));
  CALL(ex_get_ids(file, EX_SIDE_SET, side_set_ids.data()));
  Write<LO> side_class_ids_w(mesh->nents(dim - 1), -1);
  auto sides_are_exposed = mark_exposed_sides(mesh);
  classify_sides_by_exposure(mesh, sides_are_exposed);
  Write<I8> side_class_dims_w =
      deep_copy(mesh->get_array<I8>(dim - 1, "class_dim"));
  auto exposed_sides2side = collect_marked(sides_are_exposed);
  map_value_into(0, exposed_sides2side, side_class_ids_w);
  if ((classify_with & NODE_SETS) && init_params.num_node_sets) {
    int max_side_set_id = 0;
    if ((classify_with & SIDE_SETS) && side_set_ids.size()) {
      max_side_set_id =
          *std::max_element(side_set_ids.begin(), side_set_ids.end());
    }
    std::vector<int> node_set_ids(std::size_t(init_params.num_node_sets));
    CALL(ex_get_ids(file, EX_NODE_SET, node_set_ids.data()));
    std::vector<char> names_memory;
    std::vector<char*> name_ptrs;
    setup_names(int(init_params.num_node_sets), names_memory, name_ptrs);
    CALL(ex_get_names(file, EX_NODE_SET, name_ptrs.data()));
    for (size_t i = 0; i < node_set_ids.size(); ++i) {
      int nentries, ndist_factors;
      CALL(ex_get_set_param(
          file, EX_NODE_SET, node_set_ids[i], &nentries, &ndist_factors));
      if (verbose) {
        std::cout << "P" << mesh->comm()->rank() << ": node set " << node_set_ids[i] << " has " << nentries
                  << " nodes" << std::endl;
      }
      HostWrite<LO> h_set_nodes2nodes(nentries);
      CALL(ex_get_set(file, EX_NODE_SET, node_set_ids[i],
          h_set_nodes2nodes.data(), nullptr));
      auto set_nodes2nodes =
          subtract_from_each(LOs(h_set_nodes2nodes.write()), 1);
      auto nodes_are_in_set = mark_image(set_nodes2nodes, mesh->nverts());
      auto sides_are_in_set =
          mark_up_all(mesh, VERT, dim - 1, nodes_are_in_set);
      auto set_sides2side = collect_marked(sides_are_in_set);
      auto surface_id = node_set_ids[i] + max_side_set_id;
      if (verbose) {
        std::cout << "P" << mesh->comm()->rank() << ": node set #" << node_set_ids[i] << " \"" << name_ptrs[i]
                  << "\" will be surface " << surface_id << std::endl;
      }
      map_value_into(surface_id, set_sides2side, side_class_ids_w);
      map_value_into(I8(dim - 1), set_sides2side, side_class_dims_w);
      mesh->class_sets[name_ptrs[i]].push_back({I8(dim - 1), surface_id});
    }
  }
  if (classify_with & SIDE_SETS) {
    std::vector<char> names_memory;
    std::vector<char*> name_ptrs;
    setup_names(int(init_params.num_side_sets), names_memory, name_ptrs);
    CALL(ex_get_names(file, EX_SIDE_SET, name_ptrs.data()));
    for (size_t i = 0; i < side_set_ids.size(); ++i) {
      int nentries, ndist_factors;
      CALL(ex_get_set_param(
          file, EX_SIDE_SET, side_set_ids[i], &nentries, &ndist_factors));
      if (verbose) {
        std::cout << "P" << mesh->comm()->rank() << ": side set #" << side_set_ids[i] << " \"" << name_ptrs[i]
                  << "\" has " << nentries << " sides, will be surface "
                  << side_set_ids[i] << std::endl;
      }
      HostWrite<LO> h_set_sides2elem(nentries);
      HostWrite<LO> h_set_sides2local(nentries);
      CALL(ex_get_set(file, EX_SIDE_SET, side_set_ids[i],
          h_set_sides2elem.data(), h_set_sides2local.data()));
      auto set_sides2elem =
          subtract_from_each(LOs(h_set_sides2elem.write()), 1);
      auto set_sides2local = LOs(h_set_sides2local.write());
      auto elems2sides = mesh->ask_down(dim, dim - 1).ab2b;
      auto nsides_per_elem = element_degree(family, dim, dim - 1);
      auto set_sides2side_w = Write<LO>(nentries);
      auto f2 = OMEGA_H_LAMBDA(LO set_side) {
        auto elem = set_sides2elem[set_side];
        auto side_of_element =
            side_exo2osh(family, dim, set_sides2local[set_side]);
        OMEGA_H_CHECK(side_of_element != -1);
        auto side = elems2sides[elem * nsides_per_elem + side_of_element];
        set_sides2side_w[set_side] = side;
      };
      parallel_for(nentries, f2, "set_sides2side");
      auto set_sides2side = LOs(set_sides2side_w);
      auto surface_id = side_set_ids[i];
      map_value_into(surface_id, set_sides2side, side_class_ids_w);
      map_value_into(I8(dim - 1), set_sides2side, side_class_dims_w);
      mesh->class_sets[name_ptrs[i]].push_back({I8(dim - 1), surface_id});
    }
  }
  auto elem_class_ids = LOs(elem_class_ids_w);
  auto side_class_ids = LOs(side_class_ids_w);
  auto side_class_dims = Read<I8>(side_class_dims_w);
  mesh->add_tag(dim, "class_id", 1, elem_class_ids);
  mesh->add_tag(dim - 1, "class_id", 1, side_class_ids);
  mesh->set_tag(dim - 1, "class_dim", side_class_dims);
  finalize_classification(mesh);
  end_code();
}

#if defined(OMEGA_H_USE_MPI) && defined(PARALLEL_AWARE_EXODUS)
static void read_sliced_nodal_fields(Mesh* mesh, int file, int time_step,
    bool verbose, Dist slice_verts2verts, GO nodes_begin, LO nslice_nodes) {
  int num_nodal_vars;
  CALL(ex_get_variable_param(file, EX_NODAL, &num_nodal_vars));
  if (verbose) std::cout << "P" << mesh->comm()->rank() << ": " << num_nodal_vars << " nodal variables\n";
  std::vector<char> names_memory;
  std::vector<char*> name_ptrs;
  setup_names(num_nodal_vars, names_memory, name_ptrs);
  CALL(ex_get_variable_names(file, EX_NODAL, num_nodal_vars, name_ptrs.data()));
  for (int i = 0; i < num_nodal_vars; ++i) {
    auto name = name_ptrs[std::size_t(i)];
    if (verbose) std::cout << "P" << mesh->comm()->rank() << ": Loading nodal variable \"" << name << "\"\n";
    HostWrite<double> host_write(nslice_nodes);
    CALL(ex_get_partial_var(file, time_step + 1, EX_NODAL, i + 1, /*obj_id*/ 0,
        nodes_begin + 1, nslice_nodes, host_write.data()));
    auto device_write = host_write.write();
    auto slice_data = Reals(device_write);
    auto data = slice_verts2verts.exch(slice_data, 1);
    mesh->add_tag(VERT, name, 1, data);
  }
}

Mesh read_sliced(filesystem::path const& path, CommPtr comm, bool verbose, int,
    int time_step) {
  ScopedTimer timer("exodus::read");
  verbose = verbose && (comm->rank() == 0);
  auto comm_mpi = comm->get_impl();
  auto comp_ws = int(sizeof(Real));
  int io_ws = 0;
  float version;
  auto mode = EX_READ | EX_BULK_INT64_API | EX_MPIIO;
  auto file = ex_open_par(
      path.c_str(), mode, &comp_ws, &io_ws, &version, comm_mpi, MPI_INFO_NULL);
  if (file < 0)
    Omega_h_fail("can't open sliced Exodus file %s\n", path.c_str());
  ex_init_params init_params;
  CALL(ex_get_init_ext(file, &init_params));
  if (verbose) {
    std::ostringstream oss;
    oss << "P" << comm->rank() << ": init params for " << path << ":\n";
    oss << " ExodusII " << version << '\n';
    oss << " Exodus ID " << file << '\n';
    oss << " comp_ws " << comp_ws << '\n';
    oss << " io_ws " << io_ws << '\n';
    oss << " Title " << init_params.title << '\n';
    oss << " num_dim " << init_params.num_dim << '\n';
    oss << " num_nodes " << init_params.num_nodes << '\n';
    oss << " num_elem " << init_params.num_elem << '\n';
    oss << " num_elem_blk " << init_params.num_elem_blk << '\n';
    oss << " num_node_sets " << init_params.num_node_sets << '\n';
    oss << " num_side_sets " << init_params.num_side_sets;
    std::cout << oss.str() << std::endl;
  }
  auto dim = int(init_params.num_dim);
  GO nodes_begin, nodes_end;
  suggest_slices(init_params.num_nodes, comm->size(), comm->rank(),
      &nodes_begin, &nodes_end);
  auto nslice_nodes = LO(nodes_end - nodes_begin);
  HostWrite<Real> h_coord_blk[3];
  for (Int i = 0; i < dim; ++i) {
    h_coord_blk[i] = HostWrite<Real>(nslice_nodes);
  }
  nc_set_log_level(5);
  CALL(ex_get_partial_coord(file, nodes_begin + 1, nslice_nodes,
      h_coord_blk[0].data(), h_coord_blk[1].data(), h_coord_blk[2].data()));
  HostWrite<Real> h_coords(nslice_nodes * dim);
  for (LO i = 0; i < nslice_nodes; ++i) {
    for (Int j = 0; j < dim; ++j) {
      h_coords[i * dim + j] = h_coord_blk[j][i];
    }
  }
  auto slice_coords = Reals(h_coords.write());
  std::vector<int> block_ids(std::size_t(init_params.num_elem_blk));
  CALL(ex_get_ids(file, EX_ELEM_BLOCK, block_ids.data()));
  GO elems_begin, elems_end;
  suggest_slices(init_params.num_elem, comm->size(), comm->rank(), &elems_begin,
      &elems_end);
  auto nslice_elems = LO(elems_end - elems_begin);
  HostWrite<GO> h_conn;
  Write<ClassId> elem_class_ids_w(nslice_elems);
  GO total_elem_offset = 0;
  LO slice_elem_offset = 0;
  int family_int = -1;
  for (size_t i = 0; i < block_ids.size(); ++i) {
    char elem_type[MAX_STR_LENGTH + 1];
    GO nentries;
    GO nnodes_per_entry;
    GO nedges_per_entry;
    GO nfaces_per_entry;
    GO nattr_per_entry;
    CALL(ex_get_block(file, EX_ELEM_BLOCK, block_ids[i], elem_type, &nentries,
        &nnodes_per_entry, &nedges_per_entry, &nfaces_per_entry,
        &nattr_per_entry));
    if (verbose) {
      std::cout << "P" << comm->rank() << ": block " << block_ids[i] << " has " << nentries
                << " elements of type " << elem_type << '\n';
    }
    /* some pretty weird blocks from the CDFEM people... */
    if (std::string("NULL") == elem_type && nentries == 0) continue;
    int dim_from_type;
    Omega_h_Family family_from_type;
    get_elem_type_info(elem_type, &dim_from_type, &family_from_type);
    OMEGA_H_CHECK_OP(dim_from_type, ==, dim);
    if (family_int == -1) family_int = family_from_type;
    OMEGA_H_CHECK_OP(family_int, ==, family_from_type);
    auto deg = element_degree(Omega_h_Family(family_int), dim, VERT);
    OMEGA_H_CHECK_OP(nnodes_per_entry, ==, deg);
    if (!h_conn.exists())
      h_conn = decltype(h_conn)(nslice_elems * deg, "host connectivity");
    if (nedges_per_entry < 0) nedges_per_entry = 0;
    if (nfaces_per_entry < 0) nfaces_per_entry = 0;
    if (elems_end <= total_elem_offset) continue;
    if (total_elem_offset + nentries <= elems_begin) continue;
    auto block_begin = std::max(elems_begin - total_elem_offset, GO(0));
    auto block_end = std::min(elems_end - total_elem_offset, nentries);
    auto nfrom_block = LO(block_end - block_begin);
    std::vector<int> edge_conn(std::size_t(nfrom_block * nedges_per_entry));
    std::vector<int> face_conn(std::size_t(nfrom_block * nfaces_per_entry));
    CALL(ex_get_partial_conn(file, EX_ELEM_BLOCK, block_ids[i], block_begin + 1,
        nfrom_block, h_conn.data() + slice_elem_offset * nnodes_per_entry,
        edge_conn.data(), face_conn.data()));
    auto region_id = block_ids[i];
    auto f0 = OMEGA_H_LAMBDA(LO entry) {
      elem_class_ids_w[slice_elem_offset + entry] = region_id;
    };
    parallel_for(nfrom_block, f0, "set_elem_class_ids");
    total_elem_offset += nentries;
    slice_elem_offset += nfrom_block;
  }
  OMEGA_H_CHECK_OP(total_elem_offset, ==, init_params.num_elem);
  OMEGA_H_CHECK_OP(slice_elem_offset, ==, nslice_elems);
  Omega_h_Family family = Omega_h_Family(family_int);
  auto slice_conn = subtract_from_each(GOs(h_conn.write()), GO(1));

  Dist slice_elems2elems;
  Dist slice_verts2verts;
  LOs conn;
  assemble_slices(comm, family, dim, init_params.num_elem, elems_begin,
      slice_conn, init_params.num_nodes, nodes_begin, slice_coords,
      &slice_elems2elems, &conn, &slice_verts2verts);

  auto slice_node_globals =
      GOs{nslice_nodes, nodes_begin, 1, "slice node globals"};
  auto node_globals = slice_verts2verts.exch(slice_node_globals, 1);

  Mesh mesh(comm->library());
  build_from_elems2verts(&mesh, comm, family, dim, conn, node_globals);

  auto coords = slice_verts2verts.exch(slice_coords, dim);
  mesh.add_tag(VERT, "coordinates", dim, coords);

  classify_elements(&mesh);
  auto sides_are_exposed = mark_exposed_sides(&mesh);
  classify_sides_by_exposure(&mesh, sides_are_exposed);

  auto num_time_steps = int(ex_inquire_int(file, EX_INQ_TIME));
  if (verbose) std::cout << "P" << comm->rank() << ": " << num_time_steps << " time steps\n";
  if (num_time_steps > 0) {
    if (time_step < 0) time_step = num_time_steps - 1;
    if (verbose) std::cout << "P" << comm->rank() << ": reading time step " << time_step << std::endl;
    read_sliced_nodal_fields(&mesh, file, time_step, verbose, slice_verts2verts,
        nodes_begin, nslice_nodes);
  }
  CALL(ex_close(file));
  return mesh;
}
#else
Mesh read_sliced(filesystem::path const&, CommPtr, bool, int, int) {
  Omega_h_fail(
      "Can't read Exodus file by slices, Exodus not compiled with parallel "
      "support\n");
}
#endif

void write(
    filesystem::path const& path, Mesh* mesh, bool verbose, int classify_with) {
  begin_code("exodus::write");
  auto comp_ws = int(sizeof(Real));
  auto io_ws = comp_ws;
  auto mode = EX_CLOBBER | EX_MAPS_INT64_API;
  auto file = ex_create(path.c_str(), mode, &comp_ws, &io_ws);
  if (file < 0) Omega_h_fail("can't create Exodus file %s\n", path.c_str());
  auto title = "Omega_h " OMEGA_H_SEMVER " Exodus Output";
  std::set<LO> region_set;
  auto dim = mesh->dim();
  auto elem_class_ids = mesh->get_array<ClassId>(dim, "class_id");
  auto h_elem_class_ids = HostRead<LO>(elem_class_ids);
  for (LO i = 0; i < h_elem_class_ids.size(); ++i) {
    region_set.insert(h_elem_class_ids[i]);
  }
  auto side_class_ids = mesh->get_array<ClassId>(dim - 1, "class_id");
  auto side_class_dims = mesh->get_array<I8>(dim - 1, "class_dim");
  auto h_side_class_ids = HostRead<LO>(side_class_ids);
  auto h_side_class_dims = HostRead<I8>(side_class_dims);
  std::set<LO> surface_set;
  for (LO i = 0; i < h_side_class_ids.size(); ++i) {
    if (h_side_class_dims[i] == I8(dim - 1)) {
      surface_set.insert(h_side_class_ids[i]);
    }
  }
  auto nelem_blocks = int(region_set.size());
  auto nside_sets =
      (classify_with & exodus::SIDE_SETS) ? int(surface_set.size()) : 0;
  auto nnode_sets =
      (classify_with & exodus::NODE_SETS) ? int(surface_set.size()) : 0;
  if (verbose) {
    std::ostringstream oss;
    oss << "P" << mesh->comm()->rank() << ": init params for " << path << ":\n";
    oss << " Exodus ID " << file << '\n';
    oss << " comp_ws " << comp_ws << '\n';
    oss << " io_ws " << io_ws << '\n';
    oss << " Title " << title << '\n';
    oss << " num_dim " << dim << '\n';
    oss << " num_nodes " << mesh->nverts() << '\n';
    oss << " num_elem " << mesh->nelems() << '\n';
    oss << " num_elem_blk " << nelem_blocks << '\n';
    oss << " num_node_sets " << nnode_sets << '\n';
    oss << " num_side_sets " << nside_sets;
    std::cout << oss.str() << std::endl;
  }
  CALL(ex_put_init(file, title, dim, mesh->nverts(), mesh->nelems(),
      nelem_blocks, nnode_sets, nside_sets));
  Few<Write<Real>, 3> coord_blk;
  for (Int i = 0; i < dim; ++i) coord_blk[i] = Write<Real>(mesh->nverts());
  auto coords = mesh->coords();
  auto f0 = OMEGA_H_LAMBDA(LO i) {
    for (Int j = 0; j < dim; ++j) coord_blk[j][i] = coords[i * dim + j];
  };
  parallel_for(mesh->nverts(), f0, "copy_coords");
  HostRead<Real> h_coord_blk[3];
  for (Int i = 0; i < dim; ++i) h_coord_blk[i] = HostRead<Real>(coord_blk[i]);
  CALL(ex_put_coord(file, h_coord_blk[0].data(), h_coord_blk[1].data(),
      h_coord_blk[2].data()));
  auto all_conn = mesh->ask_elem_verts();
  auto elems2file_idx = Write<LO>(mesh->nelems());
  auto elem_file_offset = LO(0);
  for (auto block_id : region_set) {
    auto type_name = (dim == 3) ? "tetra4" : "tri3";
    auto elems_in_block = each_eq_to(elem_class_ids, block_id);
    auto block_elems2elem = collect_marked(elems_in_block);
    auto nblock_elems = block_elems2elem.size();
    if (verbose) {
      std::cout << "P" << mesh->comm()->rank() <<  ": element block " << block_id << " has " << nblock_elems
                << " of type " << type_name << '\n';
    }
    auto deg = element_degree(mesh->family(), dim, VERT);
    CALL(ex_put_block(
        file, EX_ELEM_BLOCK, block_id, type_name, nblock_elems, deg, 0, 0, 0));
    auto block_conn = read(unmap(block_elems2elem, all_conn, deg));
    auto block_conn_ex = add_to_each(block_conn, 1);
    auto h_block_conn = HostRead<LO>(block_conn_ex);
    CALL(ex_put_conn(
        file, EX_ELEM_BLOCK, block_id, h_block_conn.data(), nullptr, nullptr));
    auto f = OMEGA_H_LAMBDA(LO block_elem) {
      elems2file_idx[block_elems2elem[block_elem]] =
          elem_file_offset + block_elem;
    };
    parallel_for(nblock_elems, f);
    elem_file_offset += nblock_elems;
  }
  if (classify_with) {
    for (auto set_id : surface_set) {
      auto sides_in_set = land_each(each_eq_to(side_class_ids, set_id),
          each_eq_to(side_class_dims, I8(dim - 1)));
      if (classify_with & exodus::SIDE_SETS) {
        auto set_sides2side = collect_marked(sides_in_set);
        auto nset_sides = set_sides2side.size();
        if (verbose) {
          std::cout << "P" << mesh->comm()->rank() << ": side set " << set_id << " has " << nset_sides
                    << " sides\n";
        }
        auto sides2elems = mesh->ask_up(dim - 1, dim);
        Write<int> set_sides2elem(nset_sides);
        Write<int> set_sides2local(nset_sides);
        auto f1 = OMEGA_H_LAMBDA(LO set_side) {
          auto side = set_sides2side[set_side];
          auto side_elem = sides2elems.a2ab[side];
          auto elem = sides2elems.ab2b[side_elem];
          auto elem_in_file = elems2file_idx[elem];
          auto code = sides2elems.codes[side_elem];
          auto which_down = code_which_down(code);
          set_sides2elem[set_side] = elem_in_file + 1;
          set_sides2local[set_side] = side_osh2exo(dim, which_down);
        };
        parallel_for(nset_sides, f1, "set_sides2elem");
        auto h_set_sides2elem = HostRead<int>(set_sides2elem);
        auto h_set_sides2local = HostRead<int>(set_sides2local);
        CALL(ex_put_set_param(file, EX_SIDE_SET, set_id, nset_sides, 0));
        CALL(ex_put_set(file, EX_SIDE_SET, set_id, h_set_sides2elem.data(),
            h_set_sides2local.data()));
      }
      if (classify_with & exodus::NODE_SETS) {
        auto nodes_in_set = mark_down(mesh, dim - 1, VERT, sides_in_set);
        auto set_nodes2node = collect_marked(nodes_in_set);
        auto set_nodes2node_ex = add_to_each(set_nodes2node, 1);
        auto nset_nodes = set_nodes2node.size();
        if (verbose) {
          std::cout << "P" << mesh->comm()->rank() << ": node set " << set_id << " has " << nset_nodes
                    << " nodes\n";
        }
        auto h_set_nodes2node = HostRead<LO>(set_nodes2node_ex);
        CALL(ex_put_set_param(file, EX_NODE_SET, set_id, nset_nodes, 0));
        CALL(ex_put_set(
            file, EX_NODE_SET, set_id, h_set_nodes2node.data(), nullptr));
      }
    }
    std::vector<std::string> set_names(surface_set.size());
    for (auto& pair : mesh->class_sets) {
      auto& name = pair.first;
      for (auto& cp : pair.second) {
        if (cp.dim != I8(dim - 1)) continue;
        std::size_t index = 0;
        for (auto surface_id : surface_set) {
          if (surface_id == cp.id) {
            set_names[index] = name;
            if (verbose && (classify_with & exodus::NODE_SETS)) {
              std::cout << "P" << mesh->comm()->rank() << ": node set " << surface_id << " will be called \""
                        << name << "\"\n";
            }
            if (verbose && (classify_with & exodus::SIDE_SETS)) {
              std::cout << "P" << mesh->comm()->rank() << ": side set " << surface_id << " will be called \""
                        << name << "\"\n";
            }
          }
          ++index;
        }
      }
    }
    std::vector<char*> set_name_ptrs(surface_set.size(), nullptr);
    for (std::size_t i = 0; i < set_names.size(); ++i) {
      if (set_names[i].empty()) {
        std::stringstream ss;
        ss << "surface_" << i;
        set_names[i] = ss.str();
      }
      set_name_ptrs[i] = const_cast<char*>(set_names[i].c_str());
    }
    if (classify_with & exodus::NODE_SETS) {
      CALL(ex_put_names(file, EX_NODE_SET, set_name_ptrs.data()));
    }
    if (classify_with & exodus::SIDE_SETS) {
      CALL(ex_put_names(file, EX_SIDE_SET, set_name_ptrs.data()));
    }
  }
  CALL(ex_close(file));
  end_code();
}

#undef CALL

}  // end namespace exodus

}  // end namespace Omega_h
