#ifndef OMEGA_H_FILE_HPP
#define OMEGA_H_FILE_HPP

#include <iosfwd>
#include <vector>

#include <Omega_h_config.h>
#include <Omega_h_array.hpp>
#include <Omega_h_comm.hpp>
#include <Omega_h_defines.hpp>
#include <Omega_h_filesystem.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_tag.hpp>

namespace Omega_h {

OMEGA_H_DLL Mesh read_mesh_file(filesystem::path const& path, CommPtr comm);

bool is_little_endian_cpu();

#ifdef OMEGA_H_USE_LIBMESHB
namespace meshb {
void read(Mesh* mesh, std::string const& filepath);
void write(Mesh* mesh, std::string const& filepath, int version = 2);
void read_sol(
    Mesh* mesh, std::string const& filepath, std::string const& sol_name);
void write_sol(Mesh* mesh, std::string const& filepath,
    std::string const& sol_name, int version = 2);
}  // namespace meshb
#endif

#ifdef OMEGA_H_USE_SIMMODSUITE
namespace meshsim {
Mesh read(filesystem::path const& mesh, filesystem::path const& model,
    CommPtr comm);
void call_print(LOs a);
}  // namespace meshsim
#endif

#ifdef OMEGA_H_USE_SEACASEXODUS
namespace exodus {
enum ClassifyWith {
  NODE_SETS = 0x1,
  SIDE_SETS = 0x2,
};
int open(filesystem::path const& path, bool verbose = false);
void close(int exodus_file);
int get_num_time_steps(int exodus_file);
void read_mesh(int exodus_file, Mesh* mesh, bool verbose = false,
    int classify_with = NODE_SETS | SIDE_SETS);
void read_nodal_fields(int exodus_file, Mesh* mesh, int time_step,
    std::string const& prefix = "", std::string const& postfix = "",
    bool verbose = false);
void write(filesystem::path const& path, Mesh* mesh, bool verbose = false,
    int classify_with = NODE_SETS | SIDE_SETS);
Mesh read_sliced(filesystem::path const& path, CommPtr comm,
    bool verbose = false, int classify_with = NODE_SETS | SIDE_SETS,
    int time_step = -1);
}  // namespace exodus
#endif

namespace gmsh {
Mesh read(std::istream& stream, CommPtr comm);
Mesh read(filesystem::path const& filename, CommPtr comm);
void write(std::ostream& stream, Mesh* mesh);
void write(filesystem::path const& filepath, Mesh* mesh);

#ifdef OMEGA_H_USE_GMSH

/**
 * Load a parallel MSH file specified in format version 4.1 or higher.
 *
 * \param filename path to mesh prefix. If "/path/to/square" is passed,
 * then the loader looks for parts "/path/to/square_$((RANK+1)).msh"
 * \note The caller has the responsibility to initialize Gmsh library
 * before calling this function with \c ::gmsh::initialize
 */
Mesh read_parallel(filesystem::path filename, CommPtr comm);

/**
 * Write the specified mesh in MSH format version 4.1
 *
 * \note The caller has the responsibility to initialize Gmsh library
 * before calling this function with \c ::gmsh::initialize
 */
void write_parallel(filesystem::path const& filename, Mesh& mesh);
#endif  // OMEGA_H_USE_GMSH

}  // namespace gmsh

namespace vtk {
static constexpr bool do_compress = true;
static constexpr bool dont_compress = false;
#ifdef OMEGA_H_USE_ZLIB
#define OMEGA_H_DEFAULT_COMPRESS true
#else
#define OMEGA_H_DEFAULT_COMPRESS false
#endif
TagSet get_all_vtk_tags(Mesh* mesh, Int cell_dim);
TagSet get_all_vtk_tags_mix(Mesh* mesh, Int cell_dim);
void write_vtu(std::ostream& stream, Mesh* mesh, Int cell_dim,
    TagSet const& tags, bool compress = OMEGA_H_DEFAULT_COMPRESS);
void write_vtu(filesystem::path const& filename, Mesh* mesh, Int cell_dim,
    TagSet const& tags, bool compress = OMEGA_H_DEFAULT_COMPRESS);
void write_vtu(std::string const& filename, Mesh* mesh, Int cell_dim,
    bool compress = OMEGA_H_DEFAULT_COMPRESS);
void write_vtu(std::string const& filename, Mesh* mesh,
    bool compress = OMEGA_H_DEFAULT_COMPRESS);

void write_vtu(filesystem::path const& filename, Mesh* mesh, Topo_type max_type,
    bool compress = OMEGA_H_DEFAULT_COMPRESS);

void write_parallel(filesystem::path const& path, Mesh* mesh, Int cell_dim,
    TagSet const& tags, bool compress = OMEGA_H_DEFAULT_COMPRESS);
void write_parallel(std::string const& path, Mesh* mesh, Int cell_dim,
    bool compress = OMEGA_H_DEFAULT_COMPRESS);
void write_parallel(std::string const& path, Mesh* mesh,
    bool compress = OMEGA_H_DEFAULT_COMPRESS);

void read_parallel(filesystem::path const& pvtupath, CommPtr comm, Mesh* mesh);
void read_vtu(std::istream& stream, CommPtr comm, Mesh* mesh);

class Writer {
  Mesh* mesh_;
  filesystem::path root_path_;
  Int cell_dim_;
  bool compress_;
  I64 step_;
  std::streampos pvd_pos_;

 public:
  Writer();
  Writer(Writer const&) = default;
  Writer& operator=(Writer const&) = default;
  ~Writer() = default;
  Writer(filesystem::path const& root_path, Mesh* mesh, Int cell_dim = -1,
      Real restart_time = 0.0, bool compress = OMEGA_H_DEFAULT_COMPRESS);
  void write();
  void write(Real time);
  void write(Real time, TagSet const& tags);
  void write(I64 step, Real time, TagSet const& tags);
};
class FullWriter {
  std::vector<Writer> writers_;

 public:
  FullWriter() = default;
  FullWriter(filesystem::path const& root_path, Mesh* mesh,
      Real restart_time = 0.0, bool compress = OMEGA_H_DEFAULT_COMPRESS);
  void write(Real time);
  void write();
};
}  // end namespace vtk

namespace binary {

void write(filesystem::path const& path, Mesh* mesh);
Mesh read(filesystem::path const& path, Library* lib, bool strict = false);
Mesh read(filesystem::path const& path, CommPtr comm, bool strict = false);
I32 read(filesystem::path const& path, CommPtr comm, Mesh* mesh,
    bool strict = false);
I32 read_nparts(filesystem::path const& path, CommPtr comm);
I32 read_version(filesystem::path const& path, CommPtr comm);
void read_in_comm(
    filesystem::path const& path, CommPtr comm, Mesh* mesh, I32 version);

constexpr I32 latest_version = 9;

template <typename T>
void swap_bytes(T&);

template <typename T>
Read<T> swap_bytes(Read<T> array, bool needs_swapping);

template <typename T>
void write_value(std::ostream& stream, T val, bool needs_swapping);
template <typename T>
void read_value(std::istream& stream, T& val, bool needs_swapping);
template <typename T>
void write_array(std::ostream& stream, Read<T> array, bool is_compressed,
    bool needs_swapping);
template <typename T>
void read_array(std::istream& stream, Read<T>& array, bool is_compressed,
    bool needs_swapping);

void write(std::ostream& stream, std::string const& val, bool needs_swapping);
void read(std::istream& stream, std::string& val, bool needs_swapping);

void write(std::ostream& stream, Mesh* mesh);
void read(std::istream& stream, Mesh* mesh, I32 version);

#define INST_DECL(T)                                                           \
  extern template void swap_bytes(T&);                                         \
  extern template Read<T> swap_bytes(Read<T> array, bool needs_swapping);      \
  extern template void write_value(std::ostream& stream, T val, bool);         \
  extern template void read_value(std::istream& stream, T& val, bool);         \
  extern template void write_array(                                            \
      std::ostream& stream, Read<T> array, bool, bool);                        \
  extern template void read_array(                                             \
      std::istream& stream, Read<T>& array, bool, bool);
INST_DECL(I8)
INST_DECL(I32)
INST_DECL(I64)
INST_DECL(Real)
#undef INST_DECL

// for VTK compression headers
extern template void swap_bytes(std::uint64_t&);

}  // namespace binary

void write_reals_txt(filesystem::path const& filename, Reals a, Int ncomps);
void write_reals_txt(std::ostream& stream, Reals a, Int ncomps);
Reals read_reals_txt(filesystem::path const& filename, LO n, Int ncomps);
Reals read_reals_txt(std::istream& stream, LO n, Int ncomps);

}  // namespace Omega_h

#endif
