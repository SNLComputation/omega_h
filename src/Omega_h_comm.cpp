#include "Omega_h_comm.hpp"

#include <algorithm>
#include <string>

#include "Omega_h_array_ops.hpp"
#include "Omega_h_int_scan.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_profile.hpp"

#if defined(OMEGA_H_USE_CUDA) && !defined(OMEGA_H_USE_CUDA_AWARE_MPI)
#include "Omega_h_for.hpp"
#include "Omega_h_library.hpp"
#endif

namespace Omega_h {

#ifdef OMEGA_H_USE_MPI
#define CALL(f) \
{ \
  int omega_h_mpi_error = (f); \
  OMEGA_H_CHECK(MPI_SUCCESS == omega_h_mpi_error); \
}
#endif

Comm::Comm() {
#ifdef OMEGA_H_USE_MPI
  impl_ = MPI_COMM_NULL;
#endif
  library_ = nullptr;
}

#ifdef OMEGA_H_USE_MPI
Comm::Comm(Library* library_in, MPI_Comm impl_in)
    : impl_(impl_in), library_(library_in) {}

Comm::Comm(
    Library* library_in, MPI_Comm impl_in, Read<I32> srcs, Read<I32> dsts)
    : Comm(library_in, impl_in) {
  srcs_ = srcs;
  dsts_ = dsts;
  self_src_ = find_last(srcs_, rank());
  self_dst_ = find_last(dsts_, rank());
  host_srcs_ = HostRead<I32>(srcs_);
  host_dsts_ = HostRead<I32>(dsts_);
}
#else

Comm::Comm(Library* library_in, bool is_graph, bool sends_to_self)
    : library_(library_in) {
  if (is_graph) {
    if (sends_to_self) {
      srcs_ = Read<LO>({0});
      self_src_ = self_dst_ = 0;
    } else {
      srcs_ = Read<LO>({});
      self_src_ = self_dst_ = -1;
    }
    dsts_ = srcs_;
    host_srcs_ = HostRead<I32>(srcs_);
    host_dsts_ = HostRead<I32>(dsts_);
  } else {
    OMEGA_H_CHECK(!sends_to_self);
  }
}
#endif

Comm::~Comm() {
#ifdef OMEGA_H_USE_MPI
  CALL(MPI_Comm_free(&impl_));
#endif
}

Library* Comm::library() const { return library_; }

I32 Comm::rank() const {
#ifdef OMEGA_H_USE_MPI
  I32 r;
  CALL(MPI_Comm_rank(impl_, &r));
  return r;
#else
  return 0;
#endif
}

I32 Comm::size() const {
#ifdef OMEGA_H_USE_MPI
  I32 s;
  CALL(MPI_Comm_size(impl_, &s));
  return s;
#else
  return 1;
#endif
}

CommPtr Comm::dup() const {
#ifdef OMEGA_H_USE_MPI
  MPI_Comm impl2;
  CALL(MPI_Comm_dup(impl_, &impl2));
  return CommPtr(new Comm(library_, impl2));
#else
  return CommPtr(
      new Comm(library_, srcs_.exists(), srcs_.exists() && srcs_.size() == 1));
#endif
}

CommPtr Comm::split(I32 color, I32 key) const {
#ifdef OMEGA_H_USE_MPI
  MPI_Comm impl2;
  CALL(MPI_Comm_split(impl_, color, key, &impl2));
  return CommPtr(new Comm(library_, impl2));
#else
  (void)color;
  (void)key;
  return CommPtr(new Comm(library_, false, false));
#endif
}

#ifdef OMEGA_H_USE_MPI
static std::vector<int> sources_from_destinations(
    MPI_Comm comm, HostRead<I32> destinations) {
  OMEGA_H_TIME_FUNCTION;
  {
    ScopedTimer first_barrier_timer("first barrier");
    CALL(MPI_Barrier(comm));
  }
  char ignored_send_buffer = '\0';
  char ignored_recv_buffer;
  std::vector<MPI_Request> send_requests(std::size_t(destinations.size()));
  int const tag = 377;
  for (int i = 0; i < destinations.size(); ++i) {
    CALL(MPI_Issend(&ignored_send_buffer, 1, MPI_CHAR, destinations[i], tag,
        comm, &send_requests[std::size_t(i)]));
  }
  std::vector<int> sources;
  int locally_done_flag = 0;
  MPI_Request globally_done_request;
  int globally_done_flag = 0;
  while (!globally_done_flag) {
    MPI_Status status;
    int flag;
    int peer = MPI_ANY_SOURCE;
    CALL(MPI_Iprobe(peer, tag, comm, &flag, &status));
    if (flag) {
      peer = status.MPI_SOURCE;
      CALL(MPI_Recv(&ignored_recv_buffer, 1, MPI_CHAR, peer, tag, comm,
          MPI_STATUS_IGNORE));
      sources.push_back(peer);
    }
    if (locally_done_flag) {
      CALL(MPI_Test(
          &globally_done_request, &globally_done_flag, MPI_STATUS_IGNORE));
    } else {
      CALL(MPI_Testall(destinations.size(), send_requests.data(),
          &locally_done_flag, MPI_STATUSES_IGNORE));
      if (locally_done_flag) {
        CALL(MPI_Ibarrier(comm, &globally_done_request));
      }
    }
  }
  CALL(MPI_Barrier(comm));
  std::sort(sources.begin(), sources.end());
  return sources;
}
#endif

CommPtr Comm::graph(Read<I32> dsts) const {
#ifdef OMEGA_H_USE_MPI
  HostRead<I32> h_destinations(dsts);
  auto v_sources = sources_from_destinations(impl_, h_destinations);
  HostWrite<I32> h_sources(int(v_sources.size()));
  for (int i = 0; i < h_sources.size(); ++i)
    h_sources[i] = v_sources[std::size_t(i)];
  MPI_Comm impl2;
  CALL(MPI_Comm_dup(impl_, &impl2));
  return CommPtr(new Comm(library_, impl2, h_sources.write(), dsts));
#else
  return CommPtr(new Comm(library_, true, dsts.size() == 1));
#endif
}

CommPtr Comm::graph_adjacent(Read<I32> srcs, Read<I32> dsts) const {
#ifdef OMEGA_H_USE_MPI
  MPI_Comm impl2;
  CALL(MPI_Comm_dup(impl_, &impl2));
  return CommPtr(new Comm(library_, impl2, srcs, dsts));
#else
  OMEGA_H_CHECK(srcs == dsts);
  return CommPtr(new Comm(library_, true, dsts.size() == 1));
#endif
}

CommPtr Comm::graph_inverse() const {
  return graph_adjacent(destinations(), sources());
}

Read<I32> Comm::sources() const { return srcs_; }

Read<I32> Comm::destinations() const { return dsts_; }

template <typename T>
T Comm::allreduce(T x, Omega_h_Op op) const {
#ifdef OMEGA_H_USE_MPI
  CALL(MPI_Allreduce(
      MPI_IN_PLACE, &x, 1, MpiTraits<T>::datatype(), mpi_op(op), impl_));
#else
  (void)op;
#endif
  return x;
}

bool Comm::reduce_or(bool x) const {
  I8 y = x;
  y = allreduce(y, OMEGA_H_MAX);
  return static_cast<bool>(y);
}

bool Comm::reduce_and(bool x) const {
  I8 y = x;
  y = allreduce(y, OMEGA_H_MIN);
  return static_cast<bool>(y);
}

#ifdef OMEGA_H_USE_MPI
static void mpi_add_int128(void* a, void* b, int*, MPI_Datatype*) {
  Int128* a2 = static_cast<Int128*>(a);
  Int128* b2 = static_cast<Int128*>(b);
  *b2 = *b2 + *a2;
}
#endif

Int128 Comm::add_int128(Int128 x) const {
#ifdef OMEGA_H_USE_MPI
  MPI_Op op;
  int commute = true;
  CALL(MPI_Op_create(mpi_add_int128, commute, &op));
  CALL(MPI_Allreduce(MPI_IN_PLACE, &x, sizeof(Int128), MPI_PACKED, op, impl_));
  CALL(MPI_Op_free(&op));
#endif
  return x;
}

template <typename T>
T Comm::exscan(T x, Omega_h_Op op) const {
#ifdef OMEGA_H_USE_MPI
  CALL(MPI_Exscan(
      MPI_IN_PLACE, &x, 1, MpiTraits<T>::datatype(), mpi_op(op), impl_));
  if (rank() == 0) x = 0;
  return x;
#else
  (void)op;
  (void)x;
  return 0;
#endif
}

template <typename T>
void Comm::bcast(T& x, int root_rank) const {
#ifdef OMEGA_H_USE_MPI
  CALL(MPI_Bcast(&x, 1, MpiTraits<T>::datatype(), root_rank, impl_));
#else
  (void)x;
  (void)root_rank;
#endif
}

void Comm::bcast_string(std::string& s, int root_rank) const {
#ifdef OMEGA_H_USE_MPI
  I32 len = static_cast<I32>(s.length());
  bcast(len);
  s.resize(static_cast<std::size_t>(len));
  CALL(MPI_Bcast(&s[0], len, MPI_CHAR, root_rank, impl_));
#else
  (void)s;
  (void)root_rank;
#endif
}

#ifdef OMEGA_H_USE_MPI

static int Neighbor_allgather(HostRead<I32> sources, HostRead<I32> destinations,
    const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
    int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
  OMEGA_H_TIME_FUNCTION;
  static int const tag = 42;
  int const indegree = sources.size();
  int const outdegree = destinations.size();
  int typewidth;
  CALL(MPI_Type_size(sendtype, &typewidth));
  MPI_Request* sendrecvreqs = new MPI_Request[outdegree + indegree];
  for (int i = 0; i < outdegree; ++i) {
    CALL(MPI_Isend(sendbuf, sendcount, sendtype, destinations[i], tag, comm,
        sendrecvreqs + i));
  }
  for (int i = 0; i < indegree; ++i) {
    CALL(MPI_Irecv(static_cast<char*>(recvbuf) + i * typewidth, recvcount,
        recvtype, sources[i], tag, comm, sendrecvreqs + outdegree + i));
  }
  CALL(MPI_Waitall(outdegree + indegree, sendrecvreqs, MPI_STATUSES_IGNORE));
  delete[] sendrecvreqs;
  return MPI_SUCCESS;
}

static int Neighbor_alltoall(HostRead<I32> sources, HostRead<I32> destinations,
    const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
    int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
  OMEGA_H_TIME_FUNCTION;
  static int const tag = 42;
  int const indegree = sources.size();
  int const outdegree = destinations.size();
  int typewidth;
  CALL(MPI_Type_size(sendtype, &typewidth));
  MPI_Request* sendrecvreqs = new MPI_Request[outdegree + indegree];

  for (int i = 0; i < outdegree; ++i)
    CALL(MPI_Isend(static_cast<char const*>(sendbuf) + i * typewidth, sendcount,
        sendtype, destinations[i], tag, comm, sendrecvreqs + i));
  for (int i = 0; i < indegree; ++i)
    CALL(MPI_Irecv(static_cast<char*>(recvbuf) + i * typewidth, recvcount,
        recvtype, sources[i], tag, comm, sendrecvreqs + outdegree + i));
  CALL(MPI_Waitall(outdegree + indegree, sendrecvreqs, MPI_STATUSES_IGNORE));
  delete[] sendrecvreqs;
  return MPI_SUCCESS;
}

static Future<I32>::requests_type Neighbor_ialltoallv(HostRead<I32> sources,
    HostRead<I32> destinations, int width, const void* sendbuf,
    const int sdispls[], MPI_Datatype sendtype, void* recvbuf,
    const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
    int sendbuf_size = -1, int recvbuf_size = -1) {
  OMEGA_H_TIME_FUNCTION;
  int const tag = 42;
  int const indegree = sources.size();
  int const outdegree = destinations.size();
  int typewidth;
  CALL(MPI_Type_size(sendtype, &typewidth));
  Future<I32>::requests_type sendrecvreqs(static_cast<std::size_t>(outdegree + indegree));
  for (int i = 0; i < outdegree; ++i) {
    char const* const single_sendbuf =
        static_cast<char const*>(sendbuf) + sdispls[i] * typewidth * width;
    int const single_sendcount = (sdispls[i + 1] - sdispls[i]) * width;
    if (sendbuf_size != -1) {
      OMEGA_H_CHECK(static_cast<char const*>(sendbuf) <= single_sendbuf);
      OMEGA_H_CHECK(single_sendcount > 0);
      OMEGA_H_CHECK(typewidth > 0);
      OMEGA_H_CHECK(
          (single_sendbuf + single_sendcount * typewidth) <=
          (static_cast<char const*>(sendbuf) + sendbuf_size * typewidth));
    }
    CALL(MPI_Isend(single_sendbuf, single_sendcount, sendtype, destinations[i],
        tag, comm, sendrecvreqs.data() + i));
  }
  for (int i = 0; i < indegree; ++i) {
    char* const single_recvbuf =
        static_cast<char*>(recvbuf) + rdispls[i] * typewidth * width;
    int const single_recvcount = (rdispls[i + 1] - rdispls[i]) * width;
    if (recvbuf_size != -1) {
      OMEGA_H_CHECK(static_cast<char*>(recvbuf) <= single_recvbuf);
      OMEGA_H_CHECK(single_recvcount > 0);
      OMEGA_H_CHECK(typewidth > 0);
      OMEGA_H_CHECK((single_recvbuf + single_recvcount * typewidth) <=
                    (static_cast<char*>(recvbuf) + recvbuf_size * typewidth));
    }
    CALL(MPI_Irecv(single_recvbuf, single_recvcount, recvtype, sources[i], tag,
        comm, sendrecvreqs.data() + outdegree + i));
  }
  return sendrecvreqs;
}

#endif  // end ifdef OMEGA_H_USE_MPI

template <typename T>
Read<T> Comm::allgather(T x) const {
#ifdef OMEGA_H_USE_MPI
  HostWrite<T> recvbuf(srcs_.size());
  CALL(Neighbor_allgather(host_srcs_, host_dsts_, &x, 1,
      MpiTraits<T>::datatype(), nonnull(recvbuf.data()), 1,
      MpiTraits<T>::datatype(), impl_));
  return recvbuf.write();
#else
  if (srcs_.size() == 1) return Read<T>({x});
  return Read<T>({});
#endif
}

template <typename T>
Read<T> Comm::alltoall(Read<T> x) const {
#ifdef OMEGA_H_USE_MPI
  HostWrite<T> recvbuf(srcs_.size());
  HostRead<T> sendbuf(x);
  CALL(Neighbor_alltoall(host_srcs_, host_dsts_, nonnull(sendbuf.data()), 1,
      MpiTraits<T>::datatype(), nonnull(recvbuf.data()), 1,
      MpiTraits<T>::datatype(), impl_));
  return recvbuf.write();
#else
  return x;
#endif
}

#if defined(OMEGA_H_USE_CUDA) && !defined(OMEGA_H_USE_CUDA_AWARE_MPI)

template <typename T>
Read<T> self_send_part1(LO self_dst, LO self_src, Read<T>* p_sendbuf,
    Read<LO>* p_sdispls, Read<LO>* p_rdispls, Int width, LO threshold) {
  Read<T> self_data;
  if (self_dst < 0) return self_data;
  OMEGA_H_CHECK(self_src >= 0);
  auto sendbuf = *p_sendbuf;
  auto sdispls = *p_sdispls;
  auto rdispls = *p_rdispls;
  auto begin = sdispls.get(self_dst) * width;
  auto end = sdispls.get(self_dst + 1) * width;
  auto self_size = end - begin;
  if (self_size == sendbuf.size()) {
    self_data = sendbuf;
    sendbuf = Read<T>({});
    sdispls = LOs({0, 0});
    auto recvcounts_w = deep_copy(get_degrees(rdispls));
    recvcounts_w.set(self_src, 0);
    auto recvcounts = LOs(recvcounts_w);
    rdispls = offset_scan(recvcounts);
  } else {
    if (self_size * sizeof(T) < size_t(threshold)) return self_data;
    auto self_data_w = Write<T>(self_size);
    auto other_data_w = Write<T>(sendbuf.size() - self_size);
    auto f = OMEGA_H_LAMBDA(LO i) {
      if (i < begin)
        other_data_w[i] = sendbuf[i];
      else if (i < end)
        self_data_w[i - begin] = sendbuf[i];
      else
        other_data_w[i - self_size] = sendbuf[i];
    };
    parallel_for(sendbuf.size(), f, "self_send_part1");
    self_data = self_data_w;
    sendbuf = other_data_w;
    auto sendcounts_w = deep_copy(get_degrees(sdispls));
    auto recvcounts_w = deep_copy(get_degrees(rdispls));
    sendcounts_w.set(self_dst, 0);
    recvcounts_w.set(self_src, 0);
    auto sendcounts = LOs(sendcounts_w);
    auto recvcounts = LOs(recvcounts_w);
    sdispls = offset_scan(sendcounts);
    rdispls = offset_scan(recvcounts);
  }
  *p_sendbuf = sendbuf;
  *p_sdispls = sdispls;
  *p_rdispls = rdispls;
  return self_data;
}

template <typename T>
void self_send_part2(Read<T> self_data, LO self_src, Read<T>* p_recvbuf,
    Read<LO> rdispls, Int width) {
  if (!self_data.exists()) return;
  auto recvbuf = *p_recvbuf;
  if (recvbuf.size() == 0) {
    recvbuf = self_data;
  } else {
    auto begin = rdispls.get(self_src) * width;
    auto self_size = self_data.size();
    auto end = begin + self_size;
    auto recvbuf_w = Write<T>(recvbuf.size() + self_size);
    auto f = OMEGA_H_LAMBDA(LO i) {
      if (i < begin)
        recvbuf_w[i] = recvbuf[i];
      else if (i < end)
        recvbuf_w[i] = self_data[i - begin];
      else
        recvbuf_w[i] = recvbuf[i - self_size];
    };
    parallel_for(recvbuf_w.size(), f, "self_send_part2");
    recvbuf = recvbuf_w;
  }
  *p_recvbuf = recvbuf;
}

#endif

template <typename T>
Future<T> Comm::ialltoallv(Read<T> sendbuf_dev, Read<LO> sdispls_dev,
    Read<LO> rdispls_dev, Int width) const {
  ScopedTimer timer("Comm::ialltoallv");
#ifdef OMEGA_H_USE_MPI
#if defined(OMEGA_H_USE_CUDA) && !defined(OMEGA_H_USE_CUDA_AWARE_MPI)
  auto self_data = self_send_part1(self_dst_, self_src_, &sendbuf_dev,
      &sdispls_dev, &rdispls_dev, width, library_->self_send_threshold());
#endif
  HostRead<LO> sdispls(sdispls_dev);
  HostRead<LO> rdispls(rdispls_dev);
  OMEGA_H_CHECK(sendbuf_dev.size() == sdispls.last() * width);
  int nrecvd = rdispls.last() * width;
#if defined(OMEGA_H_USE_CUDA) && !defined(OMEGA_H_USE_CUDA_AWARE_MPI)
  HostWrite<T> recvbuf(nrecvd);
  HostRead<T> sendbuf(sendbuf_dev);
  OMEGA_H_CHECK(recvbuf.size() == rdispls.last() * width);
  auto reqs = Neighbor_ialltoallv(host_srcs_, host_dsts_, width,
      nonnull(sendbuf.data()), nonnull(sdispls.data()),
      MpiTraits<T>::datatype(), nonnull(recvbuf.data()),
      nonnull(rdispls.data()), MpiTraits<T>::datatype(), impl_);
  auto callback = [=](HostWrite<T> buf) -> Read<T> {
    auto recvbuf_dev = Read<T>(buf.write());
    self_send_part2(self_data, self_src_, &recvbuf_dev, rdispls_dev, width);
    return recvbuf_dev;
  };
  return {sendbuf, recvbuf, std::move(reqs), callback};
#else   // !defined(OMEGA_H_USE_CUDA) || defined(OMEGA_H_USE_CUDA_AWARE_MPI)
  Write<T> recvbuf_dev_w(nrecvd);
  OMEGA_H_CHECK(recvbuf_dev_w.size() == rdispls.last() * width);
  auto reqs = Neighbor_ialltoallv(host_srcs_, host_dsts_, width,
      nonnull(sendbuf_dev.data()), nonnull(sdispls.data()),
      MpiTraits<T>::datatype(), nonnull(recvbuf_dev_w.data()),
      nonnull(rdispls.data()), MpiTraits<T>::datatype(), impl_,
      sendbuf_dev.size(), recvbuf_dev_w.size());
  return {sendbuf_dev, recvbuf_dev_w, std::move(reqs)};
#endif  // !defined(OMEGA_H_USE_CUDA) || defined(OMEGA_H_USE_CUDA_AWARE_MPI)
#else   // !defined(OMEGA_H_USE_MPI)
  (void)sdispls_dev;
  (void)rdispls_dev;
  (void)width;
  return Future<T>(sendbuf_dev);
#endif  // !defined(OMEGA_H_USE_MPI)
}

template <typename T>
Read<T> Comm::alltoallv(Read<T> sendbuf_dev, Read<LO> sdispls_dev,
    Read<LO> rdispls_dev, Int width) const {
  ScopedTimer timer("Comm::alltoallv");
#ifdef OMEGA_H_USE_MPI
#if defined(OMEGA_H_USE_CUDA) && !defined(OMEGA_H_USE_CUDA_AWARE_MPI)
  auto self_data = self_send_part1(self_dst_, self_src_, &sendbuf_dev,
      &sdispls_dev, &rdispls_dev, width, library_->self_send_threshold());
#endif
  HostRead<LO> sdispls(sdispls_dev);
  HostRead<LO> rdispls(rdispls_dev);
  OMEGA_H_CHECK(sendbuf_dev.size() == sdispls.last() * width);
  int nrecvd = rdispls.last() * width;
#if defined(OMEGA_H_USE_CUDA) && !defined(OMEGA_H_USE_CUDA_AWARE_MPI)
  HostWrite<T> recvbuf(nrecvd);
  HostRead<T> sendbuf(sendbuf_dev);
  OMEGA_H_CHECK(recvbuf.size() == rdispls.last() * width);
  auto reqs = Neighbor_ialltoallv(host_srcs_, host_dsts_, width,
      nonnull(sendbuf.data()), nonnull(sdispls.data()),
      MpiTraits<T>::datatype(), nonnull(recvbuf.data()),
      nonnull(rdispls.data()), MpiTraits<T>::datatype(), impl_);
  CALL(MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE));
  auto recvbuf_dev = Read<T>(recvbuf.write());
  self_send_part2(self_data, self_src_, &recvbuf_dev, rdispls_dev, width);
#else   // !defined(OMEGA_H_USE_CUDA) || defined(OMEGA_H_USE_CUDA_AWARE_MPI)
  Write<T> recvbuf_dev_w(nrecvd);
  OMEGA_H_CHECK(recvbuf_dev_w.size() == rdispls.last() * width);
  auto reqs = Neighbor_ialltoallv(host_srcs_, host_dsts_, width,
      nonnull(sendbuf_dev.data()), nonnull(sdispls.data()),
      MpiTraits<T>::datatype(), nonnull(recvbuf_dev_w.data()),
      nonnull(rdispls.data()), MpiTraits<T>::datatype(), impl_,
      sendbuf_dev.size(), recvbuf_dev_w.size());
  CALL(MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE));
  Read<T> recvbuf_dev = recvbuf_dev_w;
#endif  // !defined(OMEGA_H_USE_CUDA) || defined(OMEGA_H_USE_CUDA_AWARE_MPI)
#else   // !defined(OMEGA_H_USE_MPI)
  (void)sdispls_dev;
  (void)rdispls_dev;
  (void)width;
  auto recvbuf_dev = sendbuf_dev;
#endif  // !defined(OMEGA_H_USE_MPI)
  return recvbuf_dev;
}

void Comm::barrier() const {
#ifdef OMEGA_H_USE_MPI
  CALL(MPI_Barrier(impl_));
#endif
}

#undef CALL

#define INST(T)                                                                \
  template T Comm::allreduce(T x, Omega_h_Op op) const;                        \
  template T Comm::exscan(T x, Omega_h_Op op) const;                           \
  template void Comm::bcast(T& x, int root_rank) const;                        \
  template Read<T> Comm::allgather(T x) const;                                 \
  template Read<T> Comm::alltoall(Read<T> x) const;                            \
  template Read<T> Comm::alltoallv(                                            \
      Read<T> sendbuf, Read<LO> sdispls, Read<LO> rdispls, Int width) const;   \
  template Future<T> Comm::ialltoallv(                                       \
      Read<T> sendbuf, Read<LO> sdispls, Read<LO> rdispls, Int width) const;

INST(I8)
INST(I32)
INST(I64)
INST(Real)
#undef INST

}  // end namespace Omega_h
