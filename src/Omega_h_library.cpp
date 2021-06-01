#include <Omega_h_config.h>
#include <Omega_h_cmdline.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_malloc.hpp>
#include <Omega_h_profile.hpp>
#include <Omega_h_dbg.hpp>

#include <csignal>
#include <cstdarg>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#ifdef OMEGA_H_DBG
Omega_h::Comm *DBG_COMM = 0;
bool dbg_print_global = false;
#endif

namespace Omega_h {

char* max_memory_stacktrace = nullptr;

char const* Library::static_version() { return OMEGA_H_SEMVER; }

char const* Library::version() { return static_version(); }

char const* Library::static_commit_id() { return OMEGA_H_COMMIT; }

char const* Library::commit_id() { return static_commit_id(); }

char const* Library::static_configure_options() { return OMEGA_H_CMAKE_ARGS; }

char const* Library::configure_options() { return static_configure_options(); }

#if defined(__x86_64__) || defined(_M_X64) && (!defined(OMEGA_H_USE_CUDA))
#include <xmmintrin.h>
// Intel system
static void enable_floating_point_exceptions() {
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
#ifdef __clang__
#pragma clang diagnostic pop
#endif
}
#else
static void enable_floating_point_exceptions() {
  Omega_h_fail("FPE enabled but not supported");
}
#endif

void Library::initialize(char const* head_desc, int* argc, char*** argv
#ifdef OMEGA_H_USE_MPI
    ,
    MPI_Comm comm_mpi
#endif
) {
  std::string lib_desc = OMEGA_H_SEMVER;
  if (lib_desc != head_desc) {
    std::stringstream msg;
    msg << "omega_h description string mismatch.\n";
    msg << "header says: " << head_desc << '\n';
    msg << "library says: " << lib_desc << '\n';
    std::string msg_str = msg.str();
    Omega_h::fail("%s\n", msg_str.c_str());
  }
  OMEGA_H_CHECK(argc != nullptr);
  OMEGA_H_CHECK(argv != nullptr);
  for (int ic = 0; ic < *argc; ic++) {
    argv_.push_back((*argv)[ic]);
  }
#ifdef OMEGA_H_USE_MPI
  int mpi_is_init;
  OMEGA_H_CHECK(MPI_SUCCESS == MPI_Initialized(&mpi_is_init));
  if (!mpi_is_init) {
    OMEGA_H_CHECK(MPI_SUCCESS == MPI_Init(argc, argv));
    we_called_mpi_init = true;
  } else {
    we_called_mpi_init = false;
  }
  MPI_Comm world_dup;
  MPI_Comm_dup(comm_mpi, &world_dup);
  world_ = CommPtr(new Comm(this, world_dup));

#ifdef OMEGA_H_DBG
  DBG_COMM = world_.get();
#endif

#else
  world_ = CommPtr(new Comm(this, false, false));
  self_ = CommPtr(new Comm(this, false, false));
#endif
  Omega_h::CmdLine cmdline;
  cmdline.add_flag(
      "--osh-memory", "print amount and stacktrace of max memory use");
  cmdline.add_flag(
      "--osh-time", "print amount of time spend in certain functions");
  cmdline.add_flag(
      "--osh-time-percent", "print amount of time spend in certain functions by percentage");
  auto& osh_time_chop_flag = cmdline.add_flag(
      "--osh-time-chop", "only print functions whose percent time is greater than given value (e.g. --osh-time-chop=2)");
  osh_time_chop_flag.add_arg<double>("0.0");
  cmdline.add_flag("--osh-time-with-filename", "add file name to function name in profile output");

  cmdline.add_flag("--osh-signal", "catch signals and print a stacktrace");
  cmdline.add_flag("--osh-fpe", "enable floating-point exceptions");
  cmdline.add_flag("--osh-silent", "suppress all output");
  cmdline.add_flag("--osh-pool", "use memory pooling");
  auto& self_send_flag =
      cmdline.add_flag("--osh-self-send", "control self send threshold");
  self_send_flag.add_arg<int>("value");
  auto& mpi_ranks_flag =
      cmdline.add_flag("--osh-mpi-ranks-per-node", "mpi ranks per node (for CUDA+MPI)");
  mpi_ranks_flag.add_arg<int>("value");
  if (argc && argv) {
    OMEGA_H_CHECK(cmdline.parse(world_, argc, *argv));
  }
  bool add_filename = false;
  if (cmdline.parsed("--osh-time-with-filename")) {
    add_filename = true;
  }
  double chop = cmdline.get<double>("--osh-time-chop", "0.0");
  if (cmdline.parsed("--osh-time")) {
    Omega_h::profile::global_singleton_history =
      new Omega_h::profile::History(world_, false, chop, add_filename);
  }
  if (cmdline.parsed("--osh-time-percent")) {
    Omega_h::profile::global_singleton_history =
      new Omega_h::profile::History(world_, true, chop, add_filename);
  }
  if (cmdline.parsed("--osh-fpe")) {
    enable_floating_point_exceptions();
  }
  self_send_threshold_ = 1000 * 1000;
  if (cmdline.parsed("--osh-self-send")) {
    self_send_threshold_ = cmdline.get<int>("--osh-self-send", "value");
  }
  silent_ = cmdline.parsed("--osh-silent");
#ifdef OMEGA_H_USE_KOKKOS
  if (!Kokkos::is_initialized()) {
    OMEGA_H_CHECK(argc != nullptr);
    OMEGA_H_CHECK(argv != nullptr);
    Kokkos::initialize(*argc, *argv);
    we_called_kokkos_init = true;
  } else {
    we_called_kokkos_init = false;
  }
#endif
#if defined(OMEGA_H_USE_CUDA) && defined(OMEGA_H_USE_MPI) \
  && (!defined(OMEGA_H_USE_KOKKOS))
  if (cmdline.parsed("--osh-mpi-ranks-per-node")) {
    int rank, ndevices_per_node, my_device;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    cudaGetDeviceCount(&ndevices_per_node);
    int mpi_ranks_per_node =
      cmdline.get<int>("--osh-mpi-ranks-per-node", "value");
    int local_mpi_rank = rank % mpi_ranks_per_node;
    cudaSetDevice(local_mpi_rank);
    cudaGetDevice(&my_device);
    PCOUT("ndevices_per_node= " << ndevices_per_node << " mpi_ranks_per_node= " << mpi_ranks_per_node << " local_mpi_rank= " << local_mpi_rank << std::endl);
    OMEGA_H_CHECK_OP(mpi_ranks_per_node, ==, ndevices_per_node);
    OMEGA_H_CHECK_OP(my_device, ==, local_mpi_rank);
  }
#endif
  if (cmdline.parsed("--osh-signal")) Omega_h::protect();
#if defined(OMEGA_H_USE_CUDA) && (!defined(OMEGA_H_USE_KOKKOS))
  // trigger lazy initialization of the CUDA runtime
  // and prevent it from polluting later timings
  cudaFree(nullptr);
#endif
  if (cmdline.parsed("--osh-pool")) enable_pooling();
}

Library::Library(Library const& other)
    : world_(other.world_),
      self_(other.self_)
#ifdef OMEGA_H_USE_MPI
      ,
      we_called_mpi_init(other.we_called_mpi_init)
#endif
#ifdef OMEGA_H_USE_KOKKOS
      ,
      we_called_kokkos_init(other.we_called_kokkos_init)
#endif
{
}

Library::~Library() {
  if (Omega_h::profile::global_singleton_history) {
    double total_runtime = now() - Omega_h::profile::global_singleton_history->start_time;
    if (world_->rank() == 0) {
      // FIXME - parallelize?
      Omega_h::profile::print_top_down_and_bottom_up(
          *Omega_h::profile::global_singleton_history, total_runtime);
    }
    Omega_h::profile::print_top_sorted(
          *Omega_h::profile::global_singleton_history, total_runtime);
    delete Omega_h::profile::global_singleton_history;
    Omega_h::profile::global_singleton_history = nullptr;
  }
  // need to destroy all Comm objects prior to MPI_Finalize()
  world_ = CommPtr();
  self_ = CommPtr();
  disable_pooling();
#ifdef OMEGA_H_USE_KOKKOS
  if (we_called_kokkos_init) {
    Kokkos::finalize();
    we_called_kokkos_init = false;
  }
#endif
#ifdef OMEGA_H_USE_MPI
  if (we_called_mpi_init) {
    OMEGA_H_CHECK(MPI_SUCCESS == MPI_Finalize());
    we_called_mpi_init = false;
  }
#endif
  delete[] Omega_h::max_memory_stacktrace;
}

CommPtr Library::world() { return world_; }

CommPtr Library::self() {
#ifdef OMEGA_H_USE_MPI
  if (!self_) {
    MPI_Comm self_dup;
    MPI_Comm_dup(MPI_COMM_SELF, &self_dup);
    self_ = CommPtr(new Comm(this, self_dup));
  }
#endif
  return self_;
}

LO Library::self_send_threshold() const { return self_send_threshold_; }

}  // end namespace Omega_h
