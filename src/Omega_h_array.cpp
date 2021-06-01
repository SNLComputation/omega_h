#include "Omega_h_array.hpp"

#include <cstring>
#include <sstream>

#include "Omega_h_for.hpp"

namespace Omega_h {

/* Several C libraries including ZLib and
   OpenMPI will throw errors when input pointers
   are NULL, even if they point to arrays of size zero. */
template <typename T>
class NonNullPtr {
  static T scratch[1];

 public:
  static T* get(T* p) { return (p == nullptr) ? scratch : p; }
};
template <typename T>
T NonNullPtr<T>::scratch[1] = {0};

template <typename T>
T* nonnull(T* p) {
  return NonNullPtr<T>::get(p);
}

#ifdef OMEGA_H_USE_KOKKOS
template <typename T>
Write<T>::Write(Kokkos::View<T*> view_in) : view_(view_in) {}
#endif

template <typename T>
Write<T>::Write(LO size_in, std::string const& name_in) {
  begin_code("Write allocation");
#ifdef OMEGA_H_USE_KOKKOS
  view_ = decltype(view_)(Kokkos::ViewAllocateWithoutInitializing(name_in),
      static_cast<std::size_t>(size_in));
#else
  shared_alloc_ = decltype(shared_alloc_)(
      sizeof(T) * static_cast<std::size_t>(size_in), name_in);
#endif
  end_code();
}

template <typename T>
void fill(Write<T> a, T val) {
  auto f = OMEGA_H_LAMBDA(LO i) { a[i] = val; };
  parallel_for(a.size(), f, "Write(size,value)");
}

template <typename T>
Write<T>::Write(LO size_in, T value, std::string const& name_in)
    : Write<T>(size_in, name_in) {
  fill(*this, value);
}

template <typename T>
void fill_linear(Write<T> a, T offset, T stride) {
  auto f = OMEGA_H_LAMBDA(LO i) {
    a[i] = offset + (stride * static_cast<T>(i));
  };
  parallel_for(a.size(), f, "Write(size,offset,stride)");
}

template <typename T>
Write<T>::Write(LO size_in, T offset, T stride, std::string const& name_in)
    : Write<T>(size_in, name_in) {
  fill_linear(*this, offset, stride);
}

template <typename T>
Write<T>::Write(HostWrite<T> host_write) : Write<T>(host_write.write()) {}

template <typename T>
Write<T>::Write(std::initializer_list<T> l, std::string const& name_in)
    : Write<T>(HostWrite<T>(l, name_in)) {}

#ifdef OMEGA_H_USE_KOKKOS
template <typename T>
std::string Write<T>::name() const {
  return view_.label();
}
#else
template <typename T>
std::string const& Write<T>::name() const {
  return shared_alloc_.alloc->name;
}
#endif

template <typename T>
void Write<T>::set(LO i, T value) const {
  ScopedTimer timer("single host to device");
#ifdef OMEGA_H_CHECK_BOUNDS
    OMEGA_H_CHECK(0 <= i);
    OMEGA_H_CHECK(i < size());
#endif
#ifdef OMEGA_H_USE_CUDA
  cudaMemcpy(data() + i, &value, sizeof(T), cudaMemcpyHostToDevice);
#else
  operator[](i) = value;
#endif
}

template <typename T>
T Write<T>::get(LO i) const {
  ScopedTimer timer("single device to host");
#ifdef OMEGA_H_CHECK_BOUNDS
    OMEGA_H_CHECK(0 <= i);
    OMEGA_H_CHECK(i < size());
#endif
#ifdef OMEGA_H_USE_CUDA
  T value;
  cudaMemcpy(&value, data() + i, sizeof(T), cudaMemcpyDeviceToHost);
  return value;
#else
  return operator[](i);
#endif
}

Bytes::Bytes(Write<Byte> write) : Read<Byte>(write) {}

Bytes::Bytes(LO size_in, Byte value, std::string const& name_in)
    : Read<Byte>(size_in, value, name_in) {}

Bytes::Bytes(std::initializer_list<Byte> l, std::string const& name_in)
    : Read<Byte>(l, name_in) {}

LOs::LOs(Write<LO> write) : Read<LO>(write) {}

LOs::LOs(LO size_in, LO value, std::string const& name_in)
    : Read<LO>(size_in, value, name_in) {}

LOs::LOs(LO size_in, LO offset, LO stride, std::string const& name_in)
    : Read<LO>(size_in, offset, stride, name_in) {}

LOs::LOs(std::initializer_list<LO> l, std::string const& name_in)
    : Read<LO>(l, name_in) {}

GOs::GOs(Write<GO> write) : Read<GO>(write) {}

GOs::GOs(LO size_in, GO value, std::string const& name_in)
    : Read<GO>(size_in, value, name_in) {}

GOs::GOs(LO size_in, GO offset, GO stride, std::string const& name_in)
    : Read<GO>(size_in, offset, stride, name_in) {}

GOs::GOs(std::initializer_list<GO> l, std::string const& name_in)
    : Read<GO>(l, name_in) {}

Reals::Reals(Write<Real> write) : Read<Real>(write) {}

Reals::Reals(LO size_in, Real value, std::string const& name_in)
    : Read<Real>(size_in, value, name_in) {}

Reals::Reals(std::initializer_list<Real> l, std::string const& name_in)
    : Read<Real>(l, name_in) {}

template <typename T>
Read<T>::Read(Write<T> write) : write_(write) {}

template <typename T>
Read<T>::Read(LO size_in, T value, std::string const& name_in)
    : Read<T>(Write<T>(size_in, value, name_in)) {}

template <typename T>
Read<T>::Read(LO size_in, T offset, T stride, std::string const& name_in)
    : Read<T>(Write<T>(size_in, offset, stride, name_in)) {}

template <typename T>
Read<T>::Read(std::initializer_list<T> l, std::string const& name_in)
    : Read<T>(HostWrite<T>(l, name_in).write()) {}

#ifdef OMEGA_H_USE_KOKKOS
template <typename T>
Kokkos::View<const T*> Read<T>::view() const {
  return Kokkos::View<const T*>(write_.view());
}
#endif

template <typename T>
T Read<T>::get(LO i) const {
  return write_.get(i);
}

template <typename T>
T Read<T>::first() const {
  return get(0);
}

template <typename T>
T Read<T>::last() const {
  return get(size() - 1);
}

#ifdef OMEGA_H_USE_KOKKOS
template <class T, class... P>
inline typename Kokkos::View<T, P...>::HostMirror create_uninit_mirror(
    const Kokkos::View<T, P...>& src) {
  typedef Kokkos::View<T, P...> src_type;
  typedef typename src_type::HostMirror dst_type;
  static_assert(src_type::rank == 1, "Hardcoded for 1D Views!");
  return dst_type(Kokkos::ViewAllocateWithoutInitializing(
                      std::string("host_") + src.label()),
      src.extent(0));
}

template <class T, class... P>
inline typename Kokkos::View<T, P...>::HostMirror create_uninit_mirror_view(
    const Kokkos::View<T, P...>& src,
    typename std::enable_if<(
        std::is_same<typename Kokkos::View<T, P...>::memory_space,
            typename Kokkos::View<T, P...>::HostMirror::memory_space>::value &&
        std::is_same<typename Kokkos::View<T, P...>::data_type,
            typename Kokkos::View<T, P...>::HostMirror::data_type>::value)>::
        type* = nullptr) {
  return src;
}

template <class T, class... P>
inline typename Kokkos::View<T, P...>::HostMirror create_uninit_mirror_view(
    const Kokkos::View<T, P...>& src,
    typename std::enable_if<!(
        std::is_same<typename Kokkos::View<T, P...>::memory_space,
            typename Kokkos::View<T, P...>::HostMirror::memory_space>::value &&
        std::is_same<typename Kokkos::View<T, P...>::data_type,
            typename Kokkos::View<T, P...>::HostMirror::data_type>::value)>::
        type* = nullptr) {
  return create_uninit_mirror(src);
}
#endif

template <typename T>
HostWrite<T>::HostWrite(LO size_in, std::string const& name_in)
    : write_(size_in, name_in)
#ifdef OMEGA_H_USE_KOKKOS
      ,
      mirror_(create_uninit_mirror_view(write_.view()))
#endif
{
#if (!defined(OMEGA_H_USE_KOKKOS)) && defined(OMEGA_H_USE_CUDA)
  mirror_.reset(new T[std::size_t(write_.size())]);
#endif
}

template <typename T>
HostWrite<T>::HostWrite(
    LO size_in, T offset, T stride, std::string const& name_in)
    : HostWrite<T>(Write<T>(size_in, offset, stride, name_in)) {}

template <typename T>
HostWrite<T>::HostWrite(Write<T> write_in)
    : write_(write_in)
#ifdef OMEGA_H_USE_KOKKOS
      ,
      mirror_(create_uninit_mirror_view(write_.view()))
#endif
{
#ifdef OMEGA_H_USE_KOKKOS
  Kokkos::deep_copy(mirror_, write_.view());
#elif defined(OMEGA_H_USE_CUDA)
  mirror_.reset(new T[std::size_t(write_.size())]);
  auto const err = cudaMemcpy(mirror_.get(), write_.data(),
      std::size_t(write_.size()) * sizeof(T), cudaMemcpyDeviceToHost);
  OMEGA_H_CHECK(err == cudaSuccess);
#endif
}

template <typename T>
HostWrite<T>::HostWrite(std::initializer_list<T> l, std::string const& name_in)
    :  // an initializer_list should never have over 2 billion items...
      HostWrite<T>(static_cast<LO>(l.size()), name_in) {
  LO i = 0;
  for (auto v : l) operator[](i++) = v;
}

template <typename T>
Write<T> HostWrite<T>::write() const {
  ScopedTimer timer("array host to device");
#ifdef OMEGA_H_USE_KOKKOS
  Kokkos::deep_copy(write_.view(), mirror_);
#elif defined(OMEGA_H_USE_CUDA)
  auto const err = cudaMemcpy(write_.data(), mirror_.get(),
      std::size_t(size()) * sizeof(T), cudaMemcpyHostToDevice);
  OMEGA_H_CHECK(err == cudaSuccess);
#endif
  return write_;
}

template <typename T>
LO HostWrite<T>::size() const OMEGA_H_NOEXCEPT {
  return write_.size();
}

template <typename T>
T* HostWrite<T>::data() const {
#ifdef OMEGA_H_USE_KOKKOS
  return mirror_.data();
#elif defined(OMEGA_H_USE_CUDA)
  return mirror_.get();
#else
  return write_.data();
#endif
}

template <typename T>
void HostWrite<T>::set(LO i, T value) {
#ifdef OMEGA_H_CHECK_BOUNDS
    OMEGA_H_CHECK(0 <= i);
    OMEGA_H_CHECK(i < size());
#endif
#ifdef OMEGA_H_USE_KOKKOS
  mirror_[i] = value;
#elif defined(OMEGA_H_USE_CUDA)
  mirror_[std::size_t(i)] = value;
#else
  write_[i] = value;
#endif
}

template <typename T>
T HostWrite<T>::get(LO i) const {
#ifdef OMEGA_H_USE_KOKKOS
  return mirror_[i];
#elif defined(OMEGA_H_USE_CUDA)
  return mirror_[std::size_t(i)];
#else
  return write_[i];
#endif
}

template <typename T>
HostRead<T>::HostRead(Read<T> read) : read_(read) {
  ScopedTimer timer("array device to host");
#ifdef OMEGA_H_USE_KOKKOS
  Kokkos::View<const T*> dev_view = read.view();
  Kokkos::View<const T*, Kokkos::HostSpace> h_view =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), read.view());
  mirror_ = h_view;
#elif defined(OMEGA_H_USE_CUDA)
  mirror_.reset(new T[std::size_t(read_.size())]);
  auto const err = cudaMemcpy(mirror_.get(), read_.data(),
      std::size_t(size()) * sizeof(T), cudaMemcpyDeviceToHost);
  OMEGA_H_CHECK(err == cudaSuccess);
#endif
}

template <typename T>
LO HostRead<T>::size() const {
  return read_.size();
}

template <typename T>
T const* HostRead<T>::data() const {
#if defined(OMEGA_H_USE_KOKKOS)
  return mirror_.data();
#elif defined(OMEGA_H_USE_CUDA)
  return mirror_.get();
#else
  return read_.data();
#endif
}

template <typename T>
T HostRead<T>::get(LO i) const {
  return operator[](i);
}

template <typename T>
T HostRead<T>::last() const {
  return get(size() - 1);
}

template <class T>
void copy_into(Read<T> a, Write<T> b) {
  OMEGA_H_TIME_FUNCTION;
  OMEGA_H_CHECK(a.size() == b.size());
  auto f = OMEGA_H_LAMBDA(LO i) { b[i] = a[i]; };
  parallel_for(b.size(), f, "copy into kernel");
}

template <class T>
Write<T> deep_copy(Read<T> a, std::string const& name) {
  OMEGA_H_TIME_FUNCTION;
  auto name2 = name.empty() ? a.name() : name;
  Write<T> b(a.size(), name2);
  copy_into(a, b);
  return b;
}

#define INST(T)                                                                \
  template T* nonnull(T*);                                                     \
  template T const* nonnull(T const*);                                         \
  template class NonNullPtr<T>;                                                \
  template class Write<T>;                                                     \
  template class Read<T>;                                                      \
  template class HostWrite<T>;                                                 \
  template class HostRead<T>;                                                  \
  template void fill(Write<T> a, T val);                                       \
  template void fill_linear(Write<T> a, T, T);                                 \
  template void copy_into(Read<T> a, Write<T> b);                              \
  template Write<T> deep_copy(Read<T> a, std::string const&);

INST(I8)
INST(I32)
INST(I64)
INST(Real)
#undef INST

}  // end namespace Omega_h
