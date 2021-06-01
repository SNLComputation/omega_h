#ifndef OMEGA_H_LIBRARY_HPP
#define OMEGA_H_LIBRARY_HPP

#include <map>

#include <Omega_h_comm.hpp>

namespace Omega_h {

class Library {
 public:
  Library(Library const&);
  inline Library() : Library(nullptr, nullptr) {}
  inline Library(int* argc, char*** argv
#ifdef OMEGA_H_USE_MPI
      ,
      MPI_Comm comm_mpi = MPI_COMM_WORLD
#endif
  ) {
    initialize(OMEGA_H_SEMVER, argc, argv
#ifdef OMEGA_H_USE_MPI
        ,
        comm_mpi
#endif
    );
  }
  ~Library();
  static char const* static_version();
  char const* version();
  static char const* static_commit_id();
  char const* commit_id();
  static char const* static_configure_options();
  char const* configure_options();
  CommPtr world();
  CommPtr self();
  void add_to_timer(std::string const& name, double nsecs);
  LO self_send_threshold() const;
  LO self_send_threshold_;
  bool silent_;
  std::vector<std::string> argv_;

 private:
  void initialize(char const* head_desc, int* argc, char*** argv
#ifdef OMEGA_H_USE_MPI
      ,
      MPI_Comm comm_mpi
#endif
  );
  CommPtr world_;
  CommPtr self_;
#ifdef OMEGA_H_USE_MPI
  bool we_called_mpi_init;
#endif
#ifdef OMEGA_H_USE_KOKKOS
  bool we_called_kokkos_init;
#endif
  std::map<std::string, double> timers;
};

extern char* max_memory_stacktrace;

}  // namespace Omega_h

#endif
