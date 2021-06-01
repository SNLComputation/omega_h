#ifndef OMEGA_H_COMM_HPP
#define OMEGA_H_COMM_HPP

#include <memory>

#include <Omega_h_mpi.h>
#include <Omega_h_array.hpp>
#include <Omega_h_defines.hpp>
#include <Omega_h_fail.hpp>
#include <Omega_h_future.hpp>
#include <Omega_h_int128.hpp>
#ifndef OMEGA_H_FAIL_HPP
#error "included fail but guard not defined"
#endif
#ifndef OMEGA_H_CHECK
#error "included fail but check not defined"
#endif

namespace Omega_h {

static_assert(sizeof(int) == 4, "Omega_h::Comm assumes 32-bit int");

class Library;
class Comm;

typedef std::shared_ptr<Comm> CommPtr;

class Comm {
#ifdef OMEGA_H_USE_MPI
  MPI_Comm impl_;
#endif
  Library* library_;
  Read<I32> srcs_;
  Read<I32> dsts_;
  HostRead<I32> host_srcs_;
  HostRead<I32> host_dsts_;
  LO self_src_;
  LO self_dst_;

 public:
  Comm();
#ifdef OMEGA_H_USE_MPI
  Comm(Library* library, MPI_Comm impl);
  Comm(Library* library, MPI_Comm impl, Read<I32> srcs, Read<I32> dests);
  MPI_Comm get_impl() const { return impl_; }
#else
  Comm(Library* library, bool is_graph, bool sends_to_self);
#endif
  Comm(Comm const&) = delete;
  Comm(Comm&&) = delete;
  Comm& operator=(Comm const&) = delete;
  Comm& operator=(Comm&&) = delete;
  ~Comm();
  Library* library() const;
  I32 rank() const;
  I32 size() const;
  CommPtr dup() const;
  CommPtr split(I32 color, I32 key) const;
  CommPtr graph(Read<I32> dsts) const;
  CommPtr graph_adjacent(Read<I32> srcs, Read<I32> dsts) const;
  CommPtr graph_inverse() const;
  Read<I32> sources() const;
  Read<I32> destinations() const;
  template <typename T>
  T allreduce(T x, Omega_h_Op op) const;
  bool reduce_or(bool x) const;
  bool reduce_and(bool x) const;
  Int128 add_int128(Int128 x) const;
  template <typename T>
  T exscan(T x, Omega_h_Op op) const;
  template <typename T>
  void bcast(T& x, int root_rank=0) const;
  void bcast_string(std::string& s, int root_rank=0) const;
  template <typename T>
  Read<T> allgather(T x) const;
  template <typename T>
  Read<T> alltoall(Read<T> x) const;
  template <typename T>
  Read<T> alltoallv(
      Read<T> sendbuf, Read<LO> sdispls, Read<LO> rdispls, Int width) const;
  template <typename T>
  Future<T> ialltoallv(
      Read<T> sendbuf, Read<LO> sdispls, Read<LO> rdispls, Int width) const;
  void barrier() const;
  template<typename T>
  void send(int rank, const T& x);
  template<typename T>
  void recv(int rank, T& x);
};

#ifdef OMEGA_H_USE_MPI

#ifdef OMPI_MPI_H
/* OpenMPI defines MPI_UNWEIGHTED using (void*)
 * which causes compile errors with strict
 * compile options
 */
#define OMEGA_H_MPI_UNWEIGHTED reinterpret_cast<int*>(MPI_UNWEIGHTED)
#else
#define OMEGA_H_MPI_UNWEIGHTED MPI_UNWEIGHTED
#endif

template <class T>
struct MpiTraits;

template <>
struct MpiTraits<char> {
  static MPI_Datatype datatype() { return MPI_CHAR; }
};

template <>
struct MpiTraits<I8> {
  static MPI_Datatype datatype() { return MPI_INT8_T; }
};

template <>
struct MpiTraits<I32> {
  static MPI_Datatype datatype() { return MPI_INT32_T; }
};

template <>
struct MpiTraits<I64> {
  static MPI_Datatype datatype() { return MPI_INT64_T; }
};

template <>
struct MpiTraits<double> {
  static MPI_Datatype datatype() { return MPI_DOUBLE; }
};

template <>
struct MpiTraits<unsigned> {
  static MPI_Datatype datatype() { return MPI_UNSIGNED; }
};

template <>
struct MpiTraits<unsigned long> {
  static MPI_Datatype datatype() { return MPI_UNSIGNED_LONG; }
};

inline MPI_Op mpi_op(Omega_h_Op op) {
  switch (op) {
    case OMEGA_H_MIN:
      return MPI_MIN;
    case OMEGA_H_MAX:
      return MPI_MAX;
    case OMEGA_H_SUM:
      return MPI_SUM;
  }
  OMEGA_H_NORETURN(MPI_MIN);
}


#endif


template<typename T>
void Comm::send(int rank, const T& x) {
#ifdef OMEGA_H_USE_MPI
  size_t sz = x.size();
  int omega_h_mpi_error = MPI_Send(&sz, 1, MpiTraits<size_t>::datatype(), rank, 0, impl_);
  OMEGA_H_CHECK(MPI_SUCCESS == omega_h_mpi_error);
  omega_h_mpi_error = MPI_Send(x.data(), x.size(), MpiTraits<typename T::value_type>::datatype(), rank, 0, impl_);;
  OMEGA_H_CHECK(MPI_SUCCESS == omega_h_mpi_error);
#else
  (void)rank;
  (void)x;
#endif
}

template<typename T>
void Comm::recv(int rank, T& x) {
#ifdef OMEGA_H_USE_MPI
  size_t sz = 0;
  int omega_h_mpi_error = MPI_Recv(&sz, 1, MpiTraits<size_t>::datatype(), rank, 0, impl_, MPI_STATUS_IGNORE);
  OMEGA_H_CHECK(MPI_SUCCESS == omega_h_mpi_error);
  x.resize(sz);
  omega_h_mpi_error = MPI_Recv(x.data(), x.size(), MpiTraits<typename T::value_type>::datatype(), rank, 0, impl_, MPI_STATUS_IGNORE);
  OMEGA_H_CHECK(MPI_SUCCESS == omega_h_mpi_error);
#else
  (void)rank;
  (void)x;
#endif
}

#define OMEGA_H_EXPL_INST_DECL(T)                                              \
  extern template T Comm::allreduce(T x, Omega_h_Op op) const;                 \
  extern template T Comm::exscan(T x, Omega_h_Op op) const;                    \
  extern template void Comm::bcast(T& x, int root_rank) const;                 \
  extern template Read<T> Comm::allgather(T x) const;                          \
  extern template Read<T> Comm::alltoall(Read<T> x) const;                     \
  extern template Read<T> Comm::alltoallv(                                     \
      Read<T> sendbuf, Read<LO> sdispls, Read<LO> rdispls, Int width) const;   \
  extern template Future<T> Comm::ialltoallv(                                  \
      Read<T> sendbuf, Read<LO> sdispls, Read<LO> rdispls, Int width) const;
OMEGA_H_EXPL_INST_DECL(I8)
OMEGA_H_EXPL_INST_DECL(I32)
OMEGA_H_EXPL_INST_DECL(I64)
OMEGA_H_EXPL_INST_DECL(Real)
#undef OMEGA_H_EXPL_INST_DECL

}  // end namespace Omega_h

#endif
