#include "Omega_h_vtk.hpp"
#include "Omega_h_profile.hpp"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

#ifdef OMEGA_H_USE_ZLIB
#include <zlib.h>
#endif

#include "Omega_h_array_ops.hpp"
#include "Omega_h_base64.hpp"
#include "Omega_h_build.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_tag.hpp"
#include "Omega_h_xml_lite.hpp"

#include "Omega_h_for.hpp"

namespace Omega_h {

namespace vtk {

TagSet get_all_vtk_tags(Mesh* mesh, Int cell_dim) {
  TagSet tags;
  get_all_dim_tags(mesh, VERT, &tags);
  get_all_dim_tags(mesh, cell_dim, &tags);
  tags[VERT].insert("local");
  tags[size_t(cell_dim)].insert("local");
  if (mesh->comm()->size() > 1) {
    tags[VERT].insert("owner");
    tags[size_t(cell_dim)].insert("owner");
    if (mesh->parting() == OMEGA_H_GHOSTED) {
      tags[size_t(cell_dim)].insert("vtkGhostType");
    }
  }
  return tags;
}

TagSet get_all_vtk_tags_mix(Mesh* mesh, Int cell_dim) {
  TagSet tags;
  get_all_type_tags(mesh, VERT, Topo_type::vertex, &tags);
  if (cell_dim == 3) {
    get_all_type_tags(mesh, cell_dim, Topo_type::tetrahedron, &tags);
    get_all_type_tags(mesh, cell_dim, Topo_type::hexahedron, &tags);
    get_all_type_tags(mesh, cell_dim, Topo_type::wedge, &tags);
    get_all_type_tags(mesh, cell_dim, Topo_type::pyramid, &tags);
  }
  else if (cell_dim == 2) {
    get_all_type_tags(mesh, cell_dim, Topo_type::triangle, &tags);
    get_all_type_tags(mesh, cell_dim, Topo_type::quadrilateral, &tags);
  }
  else {
    get_all_type_tags(mesh, cell_dim, Topo_type::edge, &tags);
  }
  tags[int(Topo_type::vertex)].insert("local");
  tags[size_t(cell_dim)].insert("local");
  return tags;
}

/* start of C++ ritual dance to print a string based on
   type properties */

template <>
struct IntTraits<true, 1> {
  inline static char const* name() { return "Int8"; }
};

template <>
struct IntTraits<false, 1> {
  inline static char const* name() { return "UInt8"; }
};

template <>
struct IntTraits<true, 4> {
  inline static char const* name() { return "Int32"; }
};

template <>
struct IntTraits<true, 8> {
  inline static char const* name() { return "Int64"; }
};

template <>
struct IntTraits<false, 8> {
  inline static char const* name() { return "UInt64"; }
};

template <>
struct FloatTraits<8> {
  inline static char const* name() { return "Float64"; }
};

template <typename T, typename Enable = void>
struct Traits;

template <typename T>
struct Traits<T, typename std::enable_if<std::is_integral<T>::value>::type> {
  inline static char const* name() {
    return IntTraits<std::is_signed<T>::value, sizeof(T)>::name();
  }
};

template <typename T>
struct Traits<T,
    typename std::enable_if<std::is_floating_point<T>::value>::type> {
  inline static char const* name() { return FloatTraits<sizeof(T)>::name(); }
};

/* end of C++ ritual dance to get a string based on type properties */

template <typename T>
void describe_array(std::ostream& stream, std::string const& name, Int ncomps) {
  stream << "type=\"" << Traits<T>::name() << "\"";
  stream << " Name=\"" << name << "\"";
  stream << " NumberOfComponents=\"" << ncomps << "\"";
  stream << " format=\"binary\"";
}

static bool read_array_start_tag(std::istream& stream, Omega_h_Type* type_out,
    std::string* name_out, Int* ncomps_out) {
  auto st = xml_lite::read_tag(stream);
  if (st.elem_name != "DataArray" || st.type != xml_lite::Tag::START) {
    OMEGA_H_CHECK(st.type == xml_lite::Tag::END);
    return false;
  }
  auto type_name = st.attribs["type"];
  if (type_name == "Int8")
    *type_out = OMEGA_H_I8;
  else if (type_name == "Int32")
    *type_out = OMEGA_H_I32;
  else if (type_name == "Int64")
    *type_out = OMEGA_H_I64;
  else if (type_name == "Float64")
    *type_out = OMEGA_H_F64;
  *name_out = st.attribs["Name"];
  *ncomps_out = std::stoi(st.attribs["NumberOfComponents"]);
  OMEGA_H_CHECK(st.attribs["format"] == "binary");
  return true;
}

template <typename T_osh, typename T_vtk>
void write_array(std::ostream& stream, std::string const& name, Int ncomps,
    Read<T_osh> array, bool compress) {
  OMEGA_H_TIME_FUNCTION;
  if (!(array.exists())) {
    Omega_h_fail("vtk::write_array: \"%s\" doesn't exist\n", name.c_str());
  }
  begin_code("header");
  stream << "<DataArray ";
  describe_array<T_vtk>(stream, name, ncomps);
  stream << ">\n";
  end_code();
  HostRead<T_osh> uncompressed(array);
  std::uint64_t uncompressed_bytes =
      sizeof(T_osh) * static_cast<uint64_t>(array.size());
  std::string enc_header;
  std::string encoded;
#ifdef OMEGA_H_USE_ZLIB
  if (compress) {
    begin_code("zlib");
    uLong source_bytes = uncompressed_bytes;
    uLong dest_bytes = ::compressBound(source_bytes);
    auto compressed = new ::Bytef[dest_bytes];
    int ret = ::compress2(compressed, &dest_bytes,
        reinterpret_cast<const ::Bytef*>(nonnull(uncompressed.data())),
        source_bytes, Z_BEST_SPEED);
    end_code();
    OMEGA_H_CHECK(ret == Z_OK);
    begin_code("base64");
    encoded = base64::encode(compressed, dest_bytes);
    delete[] compressed;
    std::uint64_t header[4] = {
        1, uncompressed_bytes, uncompressed_bytes, dest_bytes};
    enc_header = base64::encode(header, sizeof(header));
    end_code();
  } else
#else
  OMEGA_H_CHECK(!compress);
#endif
  {
    begin_code("base64 bulk");
    enc_header = base64::encode(&uncompressed_bytes, sizeof(std::uint64_t));
    encoded = base64::encode(nonnull(uncompressed.data()), uncompressed_bytes);
    end_code();
  }
  begin_code("stream bulk");
  // stream << enc_header << encoded << '\n';
  // the following three lines are 30% faster than the above line
  stream.write(enc_header.data(), std::streamsize(enc_header.length()));
  stream.write(encoded.data(), std::streamsize(encoded.length()));
  stream.write("\n", 1);
  end_code();
  begin_code("footer");
  stream << "</DataArray>\n";
  end_code();
}

template <typename T>
static Read<T> read_array(
    std::istream& stream, LO size, bool needs_swapping, bool is_compressed) {
  auto enc_both = base64::read_encoded(stream);
  std::uint64_t uncompressed_bytes;
  std::string encoded;
#ifdef OMEGA_H_USE_ZLIB
  std::uint64_t compressed_bytes;
  if (is_compressed) {
    std::uint64_t header[4];
    auto nheader_chars = base64::encoded_size(sizeof(header));
    auto enc_header = enc_both.substr(0, nheader_chars);
    base64::decode(enc_header, header, sizeof(header));
    if (needs_swapping) {
      for (std::uint64_t i = 0; i < 4; ++i) {
        binary::swap_bytes(header[i]);
      }
    }
    encoded = enc_both.substr(nheader_chars);
    uncompressed_bytes = header[2];
    compressed_bytes = header[3];
  } else
#else
  OMEGA_H_CHECK(is_compressed == false);
#endif
  {
    auto nheader_chars = base64::encoded_size(sizeof(std::uint64_t));
    auto enc_header = enc_both.substr(0, nheader_chars);
    base64::decode(enc_header, &uncompressed_bytes, sizeof(uncompressed_bytes));
    if (needs_swapping) binary::swap_bytes(uncompressed_bytes);
#ifdef OMEGA_H_USE_ZLIB
    compressed_bytes = uncompressed_bytes;
#endif
    encoded = enc_both.substr(nheader_chars);
  }
  OMEGA_H_CHECK(uncompressed_bytes == std::uint64_t(size) * sizeof(T));
  HostWrite<T> uncompressed(size);
#ifdef OMEGA_H_USE_ZLIB
  if (is_compressed) {
    auto compressed = new ::Bytef[compressed_bytes];
    base64::decode(encoded, compressed, compressed_bytes);
    uLong dest_bytes = static_cast<uLong>(uncompressed_bytes);
    uLong source_bytes = static_cast<uLong>(compressed_bytes);
    ::Bytef* uncompressed_ptr =
        reinterpret_cast< ::Bytef*>(nonnull(uncompressed.data()));
    int ret =
        ::uncompress(uncompressed_ptr, &dest_bytes, compressed, source_bytes);
    if (ret != Z_OK) {
      Omega_h_fail(
          "code %d: couln't decompress \"%s\"\n", ret, encoded.c_str());
    }
    OMEGA_H_CHECK(dest_bytes == static_cast<uLong>(uncompressed_bytes));
    delete[] compressed;
  } else
#endif
  {
    base64::decode(encoded, nonnull(uncompressed.data()), uncompressed_bytes);
  }
  return binary::swap_bytes(Read<T>(uncompressed.write()), needs_swapping);
}

void write_tag(
    std::ostream& stream, TagBase const* tag, Int space_dim, bool compress) {
  OMEGA_H_TIME_FUNCTION;
  if (is<I8>(tag)) {
    write_array(
        stream, tag->name(), tag->ncomps(), as<I8>(tag)->array(), compress);
  } else if (is<I32>(tag)) {
    write_array(
        stream, tag->name(), tag->ncomps(), as<I32>(tag)->array(), compress);
  } else if (is<I64>(tag)) {
    write_array(
        stream, tag->name(), tag->ncomps(), as<I64>(tag)->array(), compress);
  } else if (is<Real>(tag)) {
    Reals array = as<Real>(tag)->array();
    if (1 < space_dim && space_dim < 3) {
      if (tag->ncomps() == space_dim) {
        // VTK / ParaView expect vector fields to have 3 components
        // regardless of whether this is a 2D mesh or not.
        // this filter adds a 3rd zero component to any
        // fields with 2 components for 2D meshes
        write_array(stream, tag->name(), 3, resize_vectors(array, space_dim, 3),
            compress);
      } else if (tag->ncomps() == symm_ncomps(space_dim)) {
        // Likewise, ParaView has component names specially set up for
        // 3D symmetric tensors
        write_array(stream, tag->name(), symm_ncomps(3),
            resize_symms(array, space_dim, 3), compress);
      } else {
        write_array(stream, tag->name(), tag->ncomps(), array, compress);
      }
    } else {
      write_array(stream, tag->name(), tag->ncomps(), array, compress);
    }
  } else {
    Omega_h_fail("unknown tag type in write_tag");
  }
}

static bool read_tag(std::istream& stream, Mesh* mesh, Int ent_dim,
    bool needs_swapping, bool is_compressed) {
  Omega_h_Type type = OMEGA_H_I8;
  std::string name;
  Int ncomps = -1;
  if (!read_array_start_tag(stream, &type, &name, &ncomps)) {
    return false;
  }
  /* tags like "global" are set by the construction mechanism,
     and it is somewhat complex to anticipate when they exist
     so we can just remove them if they are going to be reset. */
  mesh->remove_tag(ent_dim, name);
  auto size = mesh->nents(ent_dim) * ncomps;
  if (type == OMEGA_H_I8) {
    auto array = read_array<I8>(stream, size, needs_swapping, is_compressed);
    mesh->add_tag(ent_dim, name, ncomps, array, true);
  } else if (type == OMEGA_H_I32) {
    auto array = read_array<I32>(stream, size, needs_swapping, is_compressed);
    mesh->add_tag(ent_dim, name, ncomps, array, true);
  } else if (type == OMEGA_H_I64) {
    auto array = read_array<I64>(stream, size, needs_swapping, is_compressed);
    mesh->add_tag(ent_dim, name, ncomps, array, true);
  } else {
    auto array = read_array<Real>(stream, size, needs_swapping, is_compressed);
    // undo the resizes done in write_tag()
    if (1 < mesh->dim() && mesh->dim() < 3) {
      if (ncomps == 3) {
        array = resize_vectors(array, 3, mesh->dim());
        ncomps = mesh->dim();
      } else if (ncomps == symm_ncomps(3)) {
        array = resize_symms(array, 3, mesh->dim());
        ncomps = symm_ncomps(mesh->dim());
      }
    }
    mesh->add_tag(ent_dim, name, ncomps, array, true);
  }
  auto et = xml_lite::read_tag(stream);
  OMEGA_H_CHECK(et.elem_name == "DataArray");
  OMEGA_H_CHECK(et.type == xml_lite::Tag::END);
  return true;
}

template <typename T>
static Read<T> read_known_array(std::istream& stream, std::string const& name,
    LO nents, Int ncomps, bool needs_swapping, bool is_compressed) {
  auto st = xml_lite::read_tag(stream);
  OMEGA_H_CHECK(st.elem_name == "DataArray");
  OMEGA_H_CHECK(st.type == xml_lite::Tag::START);
  OMEGA_H_CHECK(st.attribs["Name"] == name);
  OMEGA_H_CHECK(st.attribs["type"] == Traits<T>::name());
  OMEGA_H_CHECK(st.attribs["NumberOfComponents"] == std::to_string(ncomps));
  auto array =
      read_array<T>(stream, nents * ncomps, needs_swapping, is_compressed);
  auto et = xml_lite::read_tag(stream);
  OMEGA_H_CHECK(et.elem_name == "DataArray");
  OMEGA_H_CHECK(et.type == xml_lite::Tag::END);
  return array;
}

enum {
  VTK_VERTEX = 1,
  VTK_POLY_VERTEX = 2,
  VTK_LINE = 3,
  VTK_POLY_LINE = 4,
  VTK_TRIANGLE = 5,
  VTK_TRIANGLE_STRIP = 6,
  VTK_POLYGON = 7,
  VTK_PIXEL = 8,
  VTK_QUAD = 9,
  VTK_TETRA = 10,
  VTK_VOXEL = 11,
  VTK_HEXAHEDRON = 12,
  VTK_WEDGE = 13,
  VTK_PYRAMID = 14
};

static constexpr I8 vtk_type(Omega_h_Family family, Int dim) {
  return (
      family == OMEGA_H_SIMPLEX
          ? (dim == 3 ? VTK_TETRA
                      : (dim == 2 ? VTK_TRIANGLE
                                  : (dim == 1 ? VTK_LINE
                                              : (dim == 0 ? VTK_VERTEX : -1))))
          : (dim == 3
                    ? VTK_HEXAHEDRON
                    : (dim == 2 ? VTK_QUAD
                                : (dim == 1 ? VTK_LINE
                                            : (dim == 0 ? VTK_VERTEX : -1)))));
}

static constexpr I8 vtk_type(Int type) {
  return (type == 7 ? VTK_PYRAMID :
         (type == 6 ? VTK_WEDGE :
         (type == 5 ? VTK_HEXAHEDRON :
         (type == 4 ? VTK_TETRA :
         (type == 3 ? VTK_QUAD :
         (type == 2 ? VTK_TRIANGLE :
         (type == 1 ? VTK_LINE :
         (type == 0 ? VTK_VERTEX : -1))))))));
}


static void read_vtkfile_vtu_start_tag(
    std::istream& stream, bool* needs_swapping_out, bool* is_compressed_out) {
  auto st = xml_lite::read_tag(stream);
  OMEGA_H_CHECK(st.elem_name == "VTKFile");
  OMEGA_H_CHECK(st.attribs["header_type"] == Traits<std::uint64_t>::name());
  auto is_little_endian = (st.attribs["byte_order"] == "LittleEndian");
  *needs_swapping_out = (is_little_endian != is_little_endian_cpu());
  auto is_compressed = (st.attribs.count("compressor") == 1);
  *is_compressed_out = is_compressed;
}

static void write_piece_start_tag(
    std::ostream& stream, Mesh const* mesh, Int cell_dim) {
  stream << "<Piece NumberOfPoints=\"" << mesh->nverts() << "\"";
  stream << " NumberOfCells=\"" << mesh->nents(cell_dim) << "\">\n";
}

static void write_piece_start_tag_mix(
    std::ostream& stream, Mesh const* mesh, Int cell_dim) {
  stream << "<Piece NumberOfPoints=\"" << mesh->nverts_mix() << "\"";
  if (cell_dim == 3) {
    stream << " NumberOfCells=\"" << mesh->nregions_mix() << "\">\n";
  }
  else if (cell_dim == 2) {
    stream << " NumberOfCells=\"" << mesh->nfaces_mix() << "\">\n";
  }
  else {
    stream << " NumberOfCells=\"" << mesh->nedges_mix() << "\">\n";
  }
}

static void read_piece_start_tag(
    std::istream& stream, LO* nverts_out, LO* ncells_out) {
  auto st = xml_lite::read_tag(stream);
  OMEGA_H_CHECK(st.elem_name == "Piece");
  *nverts_out = std::stoi(st.attribs["NumberOfPoints"]);
  *ncells_out = std::stoi(st.attribs["NumberOfCells"]);
}

static void write_connectivity(
    std::ostream& stream, Mesh* mesh, Int cell_dim, bool compress) {
  Read<I8> types(mesh->nents(cell_dim), vtk_type(mesh->family(), cell_dim));
  write_array(stream, "types", 1, types, compress);
  LOs ev2v = mesh->ask_verts_of(cell_dim);
  auto deg = element_degree(mesh->family(), cell_dim, VERT);
  /* starts off already at the end of the first entity's adjacencies,
     increments by a constant value */
  LOs ends(mesh->nents(cell_dim), deg, deg);
  write_array(stream, "connectivity", 1, ev2v, compress);
  write_array(stream, "offsets", 1, ends, compress);
}

static void write_connectivity(
    std::ostream& stream, Mesh* mesh, Int cell_dim, Topo_type max_type, bool compress) {
  if (cell_dim == 3) {
    Read<I8> types_t(mesh->nents(Topo_type::tetrahedron), vtk_type(int(Topo_type::tetrahedron)));
    Read<I8> types_h(mesh->nents(Topo_type::hexahedron), vtk_type(int(Topo_type::hexahedron)));
    Read<I8> types_w(mesh->nents(Topo_type::wedge), vtk_type(int(Topo_type::wedge)));
    Read<I8> types_p(mesh->nents(Topo_type::pyramid), vtk_type(int(Topo_type::pyramid)));
    auto types = read(concat(read(concat(read(concat(types_t, types_h)), types_w)), types_p));

    LOs tv2v = mesh->ask_verts_of(Topo_type::tetrahedron);
    LOs hv2v = mesh->ask_verts_of(Topo_type::hexahedron);
    LOs wv2v = mesh->ask_verts_of(Topo_type::wedge);
    LOs pv2v = mesh->ask_verts_of(Topo_type::pyramid);
    auto ev2v = read(concat(read(concat(read(concat(tv2v, hv2v)), wv2v)), pv2v));

    auto deg_t = element_degree(Topo_type::tetrahedron, Topo_type::vertex);
    auto deg_h = element_degree(Topo_type::hexahedron, Topo_type::vertex);
    auto deg_w = element_degree(Topo_type::wedge, Topo_type::vertex);
    auto deg_p = element_degree(Topo_type::pyramid, Topo_type::vertex);
    LOs ends_t(mesh->nents(Topo_type::tetrahedron), deg_t, deg_t);
    int lastVal = 0;
    if (ends_t.size()) {
      lastVal = ends_t.last();
    }
    LOs ends_h(mesh->nents(Topo_type::hexahedron), lastVal+deg_h, deg_h);
    if (ends_h.size()) {
      lastVal = ends_h.last();
    }
    LOs ends_w(mesh->nents(Topo_type::wedge), lastVal+deg_w, deg_w);
    if (ends_w.size()) {
      lastVal = ends_w.last();
    }
    LOs ends_p(mesh->nents(Topo_type::pyramid), lastVal+deg_p, deg_p);
    auto ends = read(concat(read(concat(read(concat(ends_t, ends_h)), ends_w)), ends_p));

    write_array(stream, "types", 1, types, compress);
    write_array(stream, "connectivity", 1, ev2v, compress);
    write_array(stream, "offsets", 1, ends, compress);
  }
  else if (cell_dim == 2) {
    Read<I8> types_tr(mesh->nents(Topo_type::triangle), vtk_type(int(Topo_type::triangle)));
    Read<I8> types_q(mesh->nents(Topo_type::quadrilateral), vtk_type(int(Topo_type::quadrilateral)));
    auto types = read(concat(types_tr, types_q));

    LOs trv2v = mesh->ask_verts_of(Topo_type::triangle);
    LOs qv2v = mesh->ask_verts_of(Topo_type::quadrilateral);
    auto ev2v = read(concat(trv2v, qv2v));

    auto deg_tr = element_degree(Topo_type::triangle, Topo_type::vertex);
    auto deg_q = element_degree(Topo_type::quadrilateral, Topo_type::vertex);
    LOs ends_tr(mesh->nents(Topo_type::triangle), deg_tr, deg_tr);
    int lastVal = 0;
    if (ends_tr.size()) {
      lastVal = ends_tr.last();
    }
    LOs ends_q(mesh->nents(Topo_type::quadrilateral), lastVal+deg_q, deg_q);
    auto ends = read(concat(ends_tr, ends_q));

    write_array(stream, "types", 1, types, compress);
    write_array(stream, "connectivity", 1, ev2v, compress);
    write_array(stream, "offsets", 1, ends, compress);
  }
  else {
    Read<I8> types(mesh->nents(max_type), vtk_type(int(max_type)));
    LOs ev2v = mesh->ask_verts_of(max_type);
    auto deg = element_degree(max_type, Topo_type::vertex);
    LOs ends(mesh->nents(max_type), deg, deg);

    write_array(stream, "types", 1, types, compress);
    write_array(stream, "connectivity", 1, ev2v, compress);
    write_array(stream, "offsets", 1, ends, compress);
  }
}

static void read_connectivity(std::istream& stream, CommPtr comm, LO ncells,
    bool needs_swapping, bool is_compressed, Omega_h_Family* family_out,
    Int* dim_out, LOs* ev2v_out) {
  auto types = read_known_array<I8>(
      stream, "types", ncells, 1, needs_swapping, is_compressed);
  Omega_h_Family family = OMEGA_H_SIMPLEX;
  Int dim = -1;
  if (types.size()) {
    auto type = types.get(0);
    if (type == VTK_LINE) {
      dim = 1;
    } else if (type == VTK_QUAD) {
      family = OMEGA_H_HYPERCUBE;
      dim = 2;
    } else if (type == VTK_HEXAHEDRON) {
      family = OMEGA_H_HYPERCUBE;
      dim = 3;
    } else if (type == VTK_TRIANGLE) {
      family = OMEGA_H_SIMPLEX;
      dim = 2;
    } else if (type == VTK_TETRA) {
      family = OMEGA_H_SIMPLEX;
      dim = 3;
    } else {
      Omega_h_fail("Unexpected VTK type %d\n", type);
    }
  }
  family = Omega_h_Family(comm->allreduce(I32(family), OMEGA_H_MAX));
  dim = comm->allreduce(dim, OMEGA_H_MAX);
  *family_out = family;
  *dim_out = dim;
  auto deg = element_degree(family, dim, VERT);
  auto ev2v = read_known_array<LO>(
      stream, "connectivity", ncells * deg, 1, needs_swapping, is_compressed);
  *ev2v_out = ev2v;
  read_known_array<LO>(
      stream, "offsets", ncells, 1, needs_swapping, is_compressed);
}

static void write_locals(
    std::ostream& stream, Mesh* mesh, Int ent_dim, bool compress) {
  write_array(
      stream, "local", 1, Read<LO>(mesh->nents(ent_dim), 0, 1), compress);
}

static void write_owners(
    std::ostream& stream, Mesh* mesh, Int ent_dim, bool compress) {
  if (mesh->comm()->size() == 1) return;
  write_array(stream, "owner", 1, mesh->ask_owners(ent_dim).ranks, compress);
}

static void write_vtk_ghost_types(
    std::ostream& stream, Mesh* mesh, Int ent_dim, bool compress) {
  if (mesh->comm()->size() == 1) return;
  const auto owned = mesh->owned(ent_dim);
  auto ghost_types = each_eq_to(owned, static_cast<I8>(0));
  write_array<I8, std::uint8_t>(
      stream, "vtkGhostType", 1, ghost_types, compress);
}

static void write_locals_and_owners(std::ostream& stream, Mesh* mesh,
    Int ent_dim, TagSet const& tags, bool compress) {
  OMEGA_H_TIME_FUNCTION;
  if (tags[size_t(ent_dim)].count("local")) {
    write_locals(stream, mesh, ent_dim, compress);
  }
  if (tags[size_t(ent_dim)].count("owner")) {
    write_owners(stream, mesh, ent_dim, compress);
  }
}

template <typename T>
void write_p_data_array(
    std::ostream& stream, std::string const& name, Int ncomps) {
  stream << "<PDataArray ";
  describe_array<T>(stream, name, ncomps);
  stream << "/>\n";
}

static void write_p_data_array2(std::ostream& stream, std::string const& name,
    Int ncomps, Int Omega_h_Type) {
  switch (Omega_h_Type) {
    case OMEGA_H_I8:
      write_p_data_array<I8>(stream, name, ncomps);
      break;
    case OMEGA_H_I32:
      write_p_data_array<I32>(stream, name, ncomps);
      break;
    case OMEGA_H_I64:
      write_p_data_array<I64>(stream, name, ncomps);
      break;
    case OMEGA_H_F64:
      write_p_data_array<Real>(stream, name, ncomps);
      break;
  }
}

void write_p_tag(std::ostream& stream, TagBase const* tag, Int space_dim) {
  if (tag->type() == OMEGA_H_REAL) {
    if (1 < space_dim && space_dim < 3) {
      if (tag->ncomps() == space_dim) {
        write_p_data_array2(stream, tag->name(), 3, OMEGA_H_REAL);
      } else if (tag->ncomps() == symm_ncomps(space_dim)) {
        write_p_data_array2(stream, tag->name(), symm_ncomps(3), OMEGA_H_REAL);
      } else {
        write_p_data_array2(stream, tag->name(), tag->ncomps(), OMEGA_H_REAL);
      }
    } else {
      write_p_data_array2(stream, tag->name(), tag->ncomps(), OMEGA_H_REAL);
    }
  } else {
    write_p_data_array2(stream, tag->name(), tag->ncomps(), tag->type());
  }
}

static filesystem::path piece_filename(
    filesystem::path const& piecepath, I32 rank) {
  auto piece_filename = piecepath;
  piece_filename += '_';
  piece_filename += std::to_string(rank);
  piece_filename += ".vtu";
  return piece_filename;
}

static filesystem::path get_rel_step_path(I64 step) {
  auto result = filesystem::path("steps");
  result /= "step_";
  result += std::to_string(step);
  return result;
}

static filesystem::path get_step_path(
    filesystem::path const& root_path, I64 step) {
  auto result = root_path;
  result /= get_rel_step_path(step);
  return result;
}

filesystem::path get_pvtu_path(filesystem::path const& step_path) {
  auto result = step_path;
  result /= "pieces.pvtu";
  return result;
}

filesystem::path get_pvd_path(filesystem::path const& root_path) {
  auto result = root_path;
  result /= "steps.pvd";
  return result;
}

static void default_dim(Mesh* mesh, Int* cell_dim) {
  if (*cell_dim == -1) *cell_dim = mesh->dim();
}

static void verify_vtk_tagset(Mesh* mesh, Int cell_dim, TagSet const& tags) {
  for (Int dim = 0; dim < 4; ++dim) {
    if (dim == 0 || dim == cell_dim) {
      for (auto& name : tags[size_t(dim)]) {
        if (!mesh->has_tag(dim, name) && name != "local" && name != "owner" &&
            name != "vtkGhostType") {
          Omega_h_fail(
              "User requested VTK output of tag %s"
              " on %s, but that tag doesn't exist on the mesh!",
              name.c_str(), dimensional_plural_name(dim));
        }
      }
    } else if (!tags[size_t(dim)].empty()) {
      Omega_h_fail(
          "User requested VTK output of tags on %s,"
          " but only vertices and %s are output!",
          dimensional_plural_name(dim), dimensional_plural_name(cell_dim));
    }
  }
}

void write_vtkfile_vtu_start_tag(std::ostream& stream, bool compress) {
  stream << "<VTKFile type=\"UnstructuredGrid\" byte_order=\"";
  if (is_little_endian_cpu())
    stream << "LittleEndian";
  else
    stream << "BigEndian";
  stream << "\" header_type=\"";
  stream << Traits<std::uint64_t>::name();
  stream << "\"";
#ifdef OMEGA_H_USE_ZLIB
  if (compress) {
    stream << " compressor=\"vtkZLibDataCompressor\"";
  }
#else
  OMEGA_H_CHECK(!compress);
#endif
  stream << ">\n";
}

void write_vtu(std::ostream& stream, Mesh* mesh, Int cell_dim,
    TagSet const& tags, bool compress) {
  OMEGA_H_TIME_FUNCTION;
  default_dim(mesh, &cell_dim);
  verify_vtk_tagset(mesh, cell_dim, tags);
  write_vtkfile_vtu_start_tag(stream, compress);
  stream << "<UnstructuredGrid>\n";
  write_piece_start_tag(stream, mesh, cell_dim);
  stream << "<Cells>\n";
  write_connectivity(stream, mesh, cell_dim, compress);
  stream << "</Cells>\n";
  stream << "<Points>\n";
  auto coords = mesh->coords();
  write_array(stream, "coordinates", 3, resize_vectors(coords, mesh->dim(), 3),
      compress);
  stream << "</Points>\n";
  stream << "<PointData>\n";
  /* globals go first so read_vtu() knows where to find them */
  if (mesh->has_tag(VERT, "global") && tags[VERT].count("global")) {
    write_tag(stream, mesh->get_tag<GO>(VERT, "global"), mesh->dim(), compress);
  }
  write_locals_and_owners(stream, mesh, VERT, tags, compress);
  for (Int i = 0; i < mesh->ntags(VERT); ++i) {
    auto tag = mesh->get_tag(VERT, i);
    if (tag->name() != "coordinates" && tag->name() != "global" &&
        tags[VERT].count(tag->name())) {
      write_tag(stream, tag, mesh->dim(), compress);
    }
  }
  stream << "</PointData>\n";
  stream << "<CellData>\n";
  /* globals go first so read_vtu() knows where to find them */
  if (mesh->has_tag(cell_dim, "global") &&
      tags[size_t(cell_dim)].count("global")) {
    write_tag(
        stream, mesh->get_tag<GO>(cell_dim, "global"), mesh->dim(), compress);
  }
  write_locals_and_owners(stream, mesh, cell_dim, tags, compress);
  if (tags[size_t(cell_dim)].count("vtkGhostType")) {
    write_vtk_ghost_types(stream, mesh, cell_dim, compress);
  }
  for (Int i = 0; i < mesh->ntags(cell_dim); ++i) {
    auto tag = mesh->get_tag(cell_dim, i);
    if (tag->name() != "global" && tags[size_t(cell_dim)].count(tag->name())) {
      write_tag(stream, tag, mesh->dim(), compress);
    }
  }
  stream << "</CellData>\n";
  stream << "</Piece>\n";
  stream << "</UnstructuredGrid>\n";
  stream << "</VTKFile>\n";
}

void write_vtu(filesystem::path const& filename, Mesh* mesh, Topo_type max_type, bool compress) {
  auto tags = get_all_vtk_tags_mix(mesh, mesh->dim());
  OMEGA_H_TIME_FUNCTION;
  std::ofstream stream(filename.c_str());
  OMEGA_H_CHECK(stream.is_open());
  auto cell_dim = mesh->ent_dim(max_type);
  default_dim(mesh, &cell_dim);
  write_vtkfile_vtu_start_tag(stream, compress);
  stream << "<UnstructuredGrid>\n";
  write_piece_start_tag_mix(stream, mesh, cell_dim);
  stream << "<Cells>\n";
  write_connectivity(stream, mesh, cell_dim, max_type, compress);
  stream << "</Cells>\n";
  stream << "<Points>\n";
  auto coords = mesh->coords_mix();

  write_array(stream, "coordinates", 3, resize_vectors(coords, mesh->dim(), 3),
      compress);
  stream << "</Points>\n";
  stream << "<PointData>\n";
  if (mesh->has_tag(VERT, "global") && tags[VERT].count("global")) {
    write_tag(stream, mesh->get_tag<GO>(VERT, "global"), mesh->dim(), compress);
  }
  for (Int i = 0; i < mesh->ntags(Topo_type::vertex); ++i) {
    auto tag = mesh->get_tag(Topo_type::vertex, i);
    if (tag->name() != "coordinates" && tag->name() != "global" &&
        tags[VERT].count(tag->name())) {
      write_tag(stream, tag, mesh->dim(), compress);
    }
  }
  stream << "</PointData>\n";
  stream << "<CellData>\n";
  if (mesh->has_tag(cell_dim, "global") &&
      tags[size_t(cell_dim)].count("global")) {
    write_tag(
        stream, mesh->get_tag<GO>(cell_dim, "global"), mesh->dim(), compress);
  }
  if (tags[size_t(cell_dim)].count("vtkGhostType")) {
    write_vtk_ghost_types(stream, mesh, cell_dim, compress);
  }
  if (cell_dim == 3) {
    for (Int i = 0; i < mesh->ntags(Topo_type::tetrahedron); ++i) {
      auto tag = mesh->get_tag(Topo_type::tetrahedron, i);
      if (tag->name() != "global" && tags[size_t(cell_dim)].count(tag->name())) {
        write_tag(stream, tag, mesh->dim(), compress);
      }
    }

    for (Int i = 0; i < mesh->ntags(Topo_type::hexahedron); ++i) {
      auto tag = mesh->get_tag(Topo_type::hexahedron, i);
      if (tag->name() != "global" && tags[size_t(cell_dim)].count(tag->name())) {
        write_tag(stream, tag, mesh->dim(), compress);
      }
    }

    for (Int i = 0; i < mesh->ntags(Topo_type::wedge); ++i) {
      auto tag = mesh->get_tag(Topo_type::wedge, i);
      if (tag->name() != "global" && tags[size_t(cell_dim)].count(tag->name())) {
        write_tag(stream, tag, mesh->dim(), compress);
      }
    }

    for (Int i = 0; i < mesh->ntags(Topo_type::pyramid); ++i) {
      auto tag = mesh->get_tag(Topo_type::pyramid, i);
      if (tag->name() != "global" && tags[size_t(cell_dim)].count(tag->name())) {
        write_tag(stream, tag, mesh->dim(), compress);
      }
    }
  }
  else if (cell_dim == 2) {
    for (Int i = 0; i < mesh->ntags(Topo_type::triangle); ++i) {
      auto tag = mesh->get_tag(Topo_type::triangle, i);
      if (tag->name() != "global" && tags[size_t(cell_dim)].count(tag->name())) {
        write_tag(stream, tag, mesh->dim(), compress);
      }
    }

    for (Int i = 0; i < mesh->ntags(Topo_type::quadrilateral); ++i) {
      auto tag = mesh->get_tag(Topo_type::quadrilateral, i);
      if (tag->name() != "global" && tags[size_t(cell_dim)].count(tag->name())) {
        write_tag(stream, tag, mesh->dim(), compress);
      }
    }
  }
  else {
    for (Int i = 0; i < mesh->ntags(cell_dim); ++i) {
      auto tag = mesh->get_tag(cell_dim, i);
      if (tag->name() != "global" && tags[size_t(cell_dim)].count(tag->name())) {
        write_tag(stream, tag, mesh->dim(), compress);
      }
    }
  }
  stream << "</CellData>\n";
  stream << "</Piece>\n";
  stream << "</UnstructuredGrid>\n";
  stream << "</VTKFile>\n";
}

void read_vtu(std::istream& stream, CommPtr comm, Mesh* mesh) {
  mesh->set_comm(comm);
  mesh->set_parting(OMEGA_H_ELEM_BASED);
  read_vtu_ents(stream, mesh);
}

void read_vtu_ents(std::istream& stream, Mesh* mesh) {
  bool needs_swapping, is_compressed;
  read_vtkfile_vtu_start_tag(stream, &needs_swapping, &is_compressed);
  auto tag1 = xml_lite::read_tag(stream);
  OMEGA_H_CHECK(tag1.elem_name == "UnstructuredGrid");
  LO nverts, ncells;
  read_piece_start_tag(stream, &nverts, &ncells);
  auto tag2 = xml_lite::read_tag(stream);
  OMEGA_H_CHECK(tag2.elem_name == "Cells");
  auto comm = mesh->comm();
  Omega_h_Family family;
  Int dim;
  LOs ev2v;
  read_connectivity(stream, comm, ncells, needs_swapping, is_compressed,
      &family, &dim, &ev2v);
  mesh->set_family(family);
  mesh->set_dim(dim);
  auto tag3 = xml_lite::read_tag(stream);
  OMEGA_H_CHECK(tag3.elem_name == "Cells");
  auto tag4 = xml_lite::read_tag(stream);
  OMEGA_H_CHECK(tag4.elem_name == "Points");
  auto coords = read_known_array<Real>(
      stream, "coordinates", nverts, 3, needs_swapping, is_compressed);
  if (dim < 3) coords = resize_vectors(coords, 3, dim);
  auto tag5 = xml_lite::read_tag(stream);
  OMEGA_H_CHECK(tag5.elem_name == "Points");
  auto tag6 = xml_lite::read_tag(stream);
  OMEGA_H_CHECK(tag6.elem_name == "PointData");
  GOs vert_globals;
  if (mesh->could_be_shared(VERT)) {
    vert_globals = read_known_array<GO>(
        stream, "global", nverts, 1, needs_swapping, is_compressed);
  } else {
    vert_globals = Read<GO>(nverts, 0, 1);
  }
  build_verts_from_globals(mesh, vert_globals);
  mesh->add_tag(VERT, "coordinates", dim, coords, true);
  while (read_tag(stream, mesh, VERT, needs_swapping, is_compressed))
    ;
  mesh->remove_tag(VERT, "local");
  mesh->remove_tag(VERT, "owner");
  auto tag7 = xml_lite::read_tag(stream);
  OMEGA_H_CHECK(tag7.elem_name == "CellData");
  GOs elem_globals;
  if (mesh->could_be_shared(dim)) {
    elem_globals = read_known_array<GO>(
        stream, "global", ncells, 1, needs_swapping, is_compressed);
  } else {
    elem_globals = Read<GO>(ncells, 0, 1);
  }
  build_ents_from_elems2verts(mesh, ev2v, vert_globals, elem_globals);
  while (read_tag(stream, mesh, dim, needs_swapping, is_compressed))
    ;
  mesh->remove_tag(dim, "local");
  mesh->remove_tag(dim, "owner");
  auto tag8 = xml_lite::read_tag(stream);
  OMEGA_H_CHECK(tag8.elem_name == "Piece");
  auto tag9 = xml_lite::read_tag(stream);
  OMEGA_H_CHECK(tag9.elem_name == "UnstructuredGrid");
  auto tag10 = xml_lite::read_tag(stream);
  OMEGA_H_CHECK(tag10.elem_name == "VTKFile");
}

void write_vtu(filesystem::path const& filename, Mesh* mesh, Int cell_dim,
    TagSet const& tags, bool compress) {
  std::ofstream file(filename.c_str());
  OMEGA_H_CHECK(file.is_open());
  ask_for_mesh_tags(mesh, tags);
  write_vtu(file, mesh, cell_dim, tags, compress);
}

void write_vtu(
    std::string const& filename, Mesh* mesh, Int cell_dim, bool compress) {
  default_dim(mesh, &cell_dim);
  write_vtu(
      filename, mesh, cell_dim, get_all_vtk_tags(mesh, cell_dim), compress);
}

void write_vtu(std::string const& filename, Mesh* mesh, bool compress) {
  write_vtu(filename, mesh, mesh->dim(), compress);
}

void write_pvtu(std::ostream& stream, Mesh* mesh, Int cell_dim,
    filesystem::path const& piecepath, TagSet const& tags) {
  OMEGA_H_TIME_FUNCTION;
  ask_for_mesh_tags(mesh, tags);
  stream << "<VTKFile type=\"PUnstructuredGrid\">\n";
  stream << "<PUnstructuredGrid";
  stream << " GhostLevel=\"";
  if (mesh->parting() == OMEGA_H_GHOSTED) {
    stream << mesh->nghost_layers();
  } else {
    if (mesh->parting() == OMEGA_H_VERT_BASED && can_print(mesh) == 0) {
      std::cerr << "WARNING: a pvtu file may be written from a "
                   "vertex-partitioned mesh, but NOT read back in\n";
    }
    stream << Int(0);
  }
  stream << "\">\n";
  stream << "<PPoints>\n";
  write_p_data_array<Real>(stream, "coordinates", 3);
  stream << "</PPoints>\n";
  stream << "<PPointData>\n";
  if (mesh->has_tag(VERT, "global") && tags[VERT].count("global")) {
    write_p_tag(stream, mesh->get_tag<GO>(VERT, "global"), mesh->dim());
  }
  if (tags[0].count("local")) {
    write_p_data_array2(stream, "local", 1, OMEGA_H_I32);
  }
  if (mesh->comm()->size() > 1 && tags[0].count("owner")) {
    write_p_data_array2(stream, "owner", 1, OMEGA_H_I32);
  }
  for (Int i = 0; i < mesh->ntags(VERT); ++i) {
    auto tag = mesh->get_tag(VERT, i);
    if (tag->name() != "coordinates" && tag->name() != "global" &&
        tags[VERT].count(tag->name())) {
      write_p_tag(stream, tag, mesh->dim());
    }
  }
  stream << "</PPointData>\n";
  stream << "<PCellData>\n";
  if (mesh->has_tag(cell_dim, "global") &&
      tags[size_t(cell_dim)].count("global")) {
    write_p_tag(stream, mesh->get_tag<GO>(cell_dim, "global"), mesh->dim());
  }
  if (tags[size_t(cell_dim)].count("local")) {
    write_p_data_array2(stream, "local", 1, OMEGA_H_I32);
  }
  if (mesh->comm()->size() > 1 && tags[size_t(cell_dim)].count("owner")) {
    write_p_data_array2(stream, "owner", 1, OMEGA_H_I32);
  }
  if (mesh->comm()->size() > 1 &&
      tags[size_t(cell_dim)].count("vtkGhostType")) {
    write_p_data_array<std::uint8_t>(stream, "vtkGhostType", 1);
  }
  for (Int i = 0; i < mesh->ntags(cell_dim); ++i) {
    auto tag = mesh->get_tag(cell_dim, i);
    if (tag->name() != "global" && tags[size_t(cell_dim)].count(tag->name())) {
      write_p_tag(stream, tag, mesh->dim());
    }
  }
  stream << "</PCellData>\n";
  for (I32 i = 0; i < mesh->comm()->size(); ++i) {
    stream << "<Piece Source=\"" << piece_filename(piecepath, i) << "\"/>\n";
  }
  stream << "</PUnstructuredGrid>\n";
  stream << "</VTKFile>\n";
}

void write_pvtu(filesystem::path const& filename, Mesh* mesh, Int cell_dim,
    filesystem::path const& piecepath, TagSet const& tags) {
  std::ofstream file(filename.c_str());
  OMEGA_H_CHECK(file.is_open());
  write_pvtu(file, mesh, cell_dim, piecepath, tags);
}

void read_pvtu(std::istream& stream, CommPtr comm, I32* npieces_out,
    std::string* vtupath_out, Int* nghost_layers_out) {
  I32 npieces = 0;
  std::string vtupath;
  Int nghost_layers = 0;
  for (std::string line; std::getline(stream, line);) {
    xml_lite::Tag tag;
    if (!xml_lite::parse_tag(line, &tag)) continue;
    if (tag.elem_name == "Piece") {
      if (npieces == comm->rank()) {
        vtupath = tag.attribs["Source"];
      }
      ++npieces;
    }
    if (tag.elem_name == "PUnstructuredGrid") {
      if (tag.attribs.find("GhostLevel") != tag.attribs.end()) {
        auto& ghostlevelstr = tag.attribs["GhostLevel"];
        nghost_layers = std::atoi(ghostlevelstr.c_str());
      }
    }
  }
  OMEGA_H_CHECK(npieces >= 1);
  OMEGA_H_CHECK(npieces <= comm->size());
  *npieces_out = npieces;
  *vtupath_out = vtupath;
  *nghost_layers_out = nghost_layers;
}

void read_pvtu(filesystem::path const& pvtupath, CommPtr comm, I32* npieces_out,
    filesystem::path* vtupath_out, Int* nghost_layers_out) {
  auto parentpath = pvtupath.parent_path();
  std::string vtupath;
  std::ifstream stream(pvtupath.c_str());
  if (!stream.is_open()) {
    Omega_h_fail("couldn't open \"%s\"\n", pvtupath.c_str());
  }
  read_pvtu(stream, comm, npieces_out, &vtupath, nghost_layers_out);
  *vtupath_out = parentpath / vtupath;
}

void write_parallel(filesystem::path const& path, Mesh* mesh, Int cell_dim,
    TagSet const& tags, bool compress) {
  ScopedTimer timer("vtk::write_parallel");
  default_dim(mesh, &cell_dim);
  ask_for_mesh_tags(mesh, tags);
  auto const rank = mesh->comm()->rank();
  if (rank == 0) {
    filesystem::create_directory(path);
  }
  mesh->comm()->barrier();
  auto const piecesdir = path / "pieces";
  if (rank == 0) {
    filesystem::create_directory(piecesdir);
  }
  mesh->comm()->barrier();
  auto const piecepath = piecesdir / "piece";
  auto const pvtuname = get_pvtu_path(path);
  if (rank == 0) {
    auto const relative_piecepath = filesystem::path("pieces") / "piece";
    write_pvtu(pvtuname, mesh, cell_dim, relative_piecepath, tags);
  }
  write_vtu(piece_filename(piecepath, rank), mesh, cell_dim, tags, compress);
}

void write_parallel(
    std::string const& path, Mesh* mesh, Int cell_dim, bool compress) {
  default_dim(mesh, &cell_dim);
  write_parallel(
      path, mesh, cell_dim, get_all_vtk_tags(mesh, cell_dim), compress);
}

void write_parallel(std::string const& path, Mesh* mesh, bool compress) {
  write_parallel(path, mesh, mesh->dim(), compress);
}

void read_parallel(filesystem::path const& pvtupath, CommPtr comm, Mesh* mesh) {
  I32 npieces;
  filesystem::path vtupath;
  Int nghost_layers;
  read_pvtu(pvtupath, comm, &npieces, &vtupath, &nghost_layers);
  bool in_subcomm = (comm->rank() < npieces);
  auto subcomm = comm->split(I32(!in_subcomm), 0);
  if (in_subcomm) {
    std::ifstream vtustream(vtupath.c_str());
    OMEGA_H_CHECK(vtustream.is_open());
    mesh->set_comm(subcomm);
    if (nghost_layers == 0) {
      mesh->set_parting(OMEGA_H_ELEM_BASED, 0, false);
    } else {
      mesh->set_parting(OMEGA_H_GHOSTED, nghost_layers, false);
    }
    read_vtu_ents(vtustream, mesh);
  }
  mesh->set_comm(comm);
}

static char const pvd_prologue[] =
    "<VTKFile type=\"Collection\" version=\"0.1\">\n<Collection>\n";
static char const pvd_epilogue[] = "</Collection>\n</VTKFile>\n";

static std::string read_existing_pvd(
    filesystem::path const& pvdpath, Real restart_time) {
  std::ifstream file(pvdpath.c_str());
  if (!file.is_open()) return pvd_prologue;
  std::string contents;
  std::string line;
  std::getline(file, contents);  // VTKFile
  contents += '\n';
  std::getline(file, line);  // Collection
  contents += line;
  contents += '\n';
  // existing file may be corrupted somehow
  if (contents != pvd_prologue) return pvd_prologue;
  while (std::getline(file, line)) {
    xml_lite::Tag tag;
    if (!xml_lite::parse_tag(line, &tag)) break;
    if (tag.elem_name != "DataSet") break;
    if (!tag.attribs.count("timestep")) break;
    auto time = std::stod(tag.attribs["timestep"]);
    if (time >= restart_time) break;
    contents += line;
    contents += '\n';
  }
  return contents;
}

std::streampos write_initial_pvd(
    filesystem::path const& root_path, Real restart_time) {
  auto const pvdpath = get_pvd_path(root_path);
  auto content = read_existing_pvd(pvdpath, restart_time);
  std::ofstream file(pvdpath.c_str());
  OMEGA_H_CHECK(file.is_open());
  file << content;
  auto pos = file.tellp();
  file << pvd_epilogue;
  return pos;
}

void update_pvd(filesystem::path const& root_path, std::streampos* pos_inout,
    I64 step, Real time) {
  auto const pvdpath = get_pvd_path(root_path);
  std::fstream file;
  file.open(pvdpath.c_str(), std::ios::out | std::ios::in);
  OMEGA_H_CHECK(file.is_open());
  file.seekp(*pos_inout);
  file << std::scientific << std::setprecision(18);
  file << "<DataSet timestep=\"" << time << "\" part=\"0\" ";
  auto const relstep = get_rel_step_path(step);
  auto const relpvtu = get_pvtu_path(relstep);
  file << "file=\"" << relpvtu << "\"/>\n";
  *pos_inout = file.tellp();
  file << "</Collection>\n";
  file << "</VTKFile>\n";
}

void read_pvd(std::istream& stream, std::vector<Real>* times_out,
    std::vector<filesystem::path>* pvtupaths_out) {
  std::vector<Real> times;
  std::vector<filesystem::path> pvtupaths;
  for (std::string line; std::getline(stream, line);) {
    xml_lite::Tag tag;
    if (!xml_lite::parse_tag(line, &tag)) continue;
    if (tag.elem_name != "DataSet") continue;
    times.push_back(std::stod(tag.attribs["timestep"]));
    pvtupaths.push_back(tag.attribs["file"]);
  }
  *times_out = times;
  *pvtupaths_out = pvtupaths;
}

void read_pvd(filesystem::path const& pvdpath, std::vector<Real>* times_out,
    std::vector<filesystem::path>* pvtupaths_out) {
  std::vector<Real> times;
  std::vector<filesystem::path> pvtupaths;
  std::ifstream pvdstream(pvdpath.c_str());
  if (!pvdstream.is_open()) {
    Omega_h_fail("Couldn't open \"%s\"\n", pvdpath.c_str());
  }
  read_pvd(pvdstream, &times, &pvtupaths);
  auto const parentpath = pvdpath.parent_path();
  for (auto& pvtupath : pvtupaths) pvtupath = parentpath / pvtupath;
  *times_out = times;
  *pvtupaths_out = pvtupaths;
}

Writer::Writer()
    : mesh_(nullptr),
      root_path_("/not-set"),
      cell_dim_(-1),
      compress_(OMEGA_H_DEFAULT_COMPRESS),
      step_(-1),
      pvd_pos_(0) {}

Writer::Writer(filesystem::path const& root_path, Mesh* mesh, Int cell_dim,
    Real restart_time, bool compress)
    : mesh_(mesh),
      root_path_(root_path),
      cell_dim_(cell_dim),
      compress_(compress),
      step_(0),
      pvd_pos_(0) {
  default_dim(mesh_, &cell_dim_);
  auto const comm = mesh->comm();
  auto const rank = comm->rank();
  if (rank == 0) {
    filesystem::create_directory(root_path_);
  }
  comm->barrier();
  auto const stepsdir = root_path_ / "steps";
  if (rank == 0) {
    filesystem::create_directory(stepsdir);
  }
  comm->barrier();
  if (rank == 0) {
    pvd_pos_ = write_initial_pvd(root_path_, restart_time);
  }
}

void Writer::write(I64 step, Real time, TagSet const& tags) {
  step_ = step;
  write_parallel(
      get_step_path(root_path_, step_), mesh_, cell_dim_, tags, compress_);
  if (mesh_->comm()->rank() == 0) {
    update_pvd(root_path_, &pvd_pos_, step_, time);
  }
}

void Writer::write(Real time, TagSet const& tags) {
  this->write(step_, time, tags);
  ++step_;
}

void Writer::write(Real time) {
  this->write(time, get_all_vtk_tags(mesh_, cell_dim_));
}

void Writer::write() { this->write(Real(step_)); }

FullWriter::FullWriter(filesystem::path const& root_path, Mesh* mesh,
    Real restart_time, bool compress) {
  auto const comm = mesh->comm();
  auto const rank = comm->rank();
  if (rank == 0) {
    filesystem::create_directory(root_path);
  }
  comm->barrier();
  for (Int i = EDGE; i <= mesh->dim(); ++i) {
    writers_.push_back(Writer(root_path / dimensional_plural_name(i), mesh, i,
        restart_time, compress));
  }
}

void FullWriter::write(Real time) {
  for (auto& writer : writers_) writer.write(time);
}

void FullWriter::write() {
  for (auto& writer : writers_) writer.write();
}

#define OMEGA_H_EXPL_INST(T)                                                   \
  template void write_p_data_array<T>(                                         \
      std::ostream & stream, std::string const& name, Int ncomps);             \
  template void write_array(std::ostream& stream, std::string const& name,     \
      Int ncomps, Read<T> array, bool compress);
OMEGA_H_EXPL_INST(I8)
OMEGA_H_EXPL_INST(I32)
OMEGA_H_EXPL_INST(I64)
OMEGA_H_EXPL_INST(Real)
#undef OMEGA_H_EXPL_INST

template void write_array<Real, std::uint8_t>(std::ostream& stream,
    std::string const& name, Int ncomps, Read<Real> array, bool compress);

}  // end namespace vtk

}  // end namespace Omega_h
