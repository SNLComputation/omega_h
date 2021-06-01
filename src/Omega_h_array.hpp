#ifndef OMEGA_H_ARRAY_HPP
#define OMEGA_H_ARRAY_HPP

#include <Omega_h_defines.hpp>
#include <Omega_h_fail.hpp>
#include <initializer_list>
#include <memory>
#ifdef OMEGA_H_USE_KOKKOS
#include <Omega_h_kokkos.hpp>
#else
#include <Omega_h_shared_alloc.hpp>
#include <string>
#endif

namespace Omega_h {

template <typename T>
T* nonnull(T* p);

template <typename T>
class HostWrite;

template <typename T>
class Write {
#ifdef OMEGA_H_USE_KOKKOS
  Kokkos::View<T*> view_;
#else
  SharedAlloc shared_alloc_;
#endif

 public:
  typedef T value_type;
  OMEGA_H_INLINE Write()
      :
#ifdef OMEGA_H_USE_KOKKOS
        view_()
#else
        shared_alloc_()
#endif
  {
  }
#ifdef OMEGA_H_USE_KOKKOS
  Write(Kokkos::View<T*> view_in);
#endif
  Write(LO size_in, std::string const& name = "");
  Write(LO size_in, T value, std::string const& name = "");
  Write(LO size_in, T offset, T stride, std::string const& name = "");
  Write(std::initializer_list<T> l, std::string const& name = "");
  Write(HostWrite<T> host_write);
  OMEGA_H_INLINE LO size() const OMEGA_H_NOEXCEPT {
#ifdef OMEGA_H_CHECK_BOUNDS
    OMEGA_H_CHECK(exists());
#endif
#ifdef OMEGA_H_USE_KOKKOS
    return static_cast<LO>(view_.size());
#else
    return static_cast<LO>(shared_alloc_.size() / sizeof(T));
#endif
  }
  OMEGA_H_DEVICE T& operator[](LO i) const OMEGA_H_NOEXCEPT {
#ifdef OMEGA_H_CHECK_BOUNDS
    OMEGA_H_CHECK_OP(0, <=, i);
    OMEGA_H_CHECK_OP(i, <, size());
#endif
#ifdef OMEGA_H_USE_KOKKOS
    return view_(i);
#else
    return data()[i];
#endif
  }
  OMEGA_H_INLINE T* data() const noexcept {
#ifdef OMEGA_H_USE_KOKKOS
    return view_.data();
#else
    return static_cast<T*>(shared_alloc_.data());
#endif
  }
#ifdef OMEGA_H_USE_KOKKOS
  OMEGA_H_INLINE Kokkos::View<T*> const& view() const { return view_; }
#endif
  void set(LO i, T value) const;
  T get(LO i) const;
#ifdef OMEGA_H_USE_KOKKOS
  OMEGA_H_INLINE long use_count() const { return view_.use_count(); }
#else
  inline int use_count() const { return shared_alloc_.alloc->use_count; }
#endif
  OMEGA_H_INLINE bool exists() const noexcept {
#if defined(OMEGA_H_USE_KOKKOS)
    return view().data() != nullptr
#if defined(KOKKOS_ENABLE_DEPRECATED_CODE) && (!defined(__CUDA_ARCH__))
           /* deprecated Kokkos behavior: zero-span views have data()==nullptr
            */
           || view().use_count() != 0
#endif
        ;
#else
    return shared_alloc_.data() != nullptr;
#endif
  }
#ifdef OMEGA_H_USE_KOKKOS
  std::string name() const;
#else
  std::string const& name() const;
#endif
  OMEGA_H_INLINE T* begin() const noexcept { return data(); }
  OMEGA_H_INLINE T* end() const OMEGA_H_NOEXCEPT { return data() + size(); }
};

template <typename T>
class Read {
  Write<T> write_;

 public:
  typedef T value_type;
  OMEGA_H_INLINE Read() {}
  Read(Write<T> write);
  Read(LO size, T value, std::string const& name = "");
  Read(LO size, T offset, T stride, std::string const& name = "");
  Read(std::initializer_list<T> l, std::string const& name = "");
  OMEGA_H_INLINE LO size() const OMEGA_H_NOEXCEPT { return write_.size(); }
  OMEGA_H_DEVICE T const& operator[](LO i) const OMEGA_H_NOEXCEPT {
#ifdef OMEGA_H_CHECK_BOUNDS
    OMEGA_H_CHECK_OP(0, <=, i);
    OMEGA_H_CHECK_OP(i, <, size());
#endif
    return write_[i];
  }
  OMEGA_H_INLINE T const* data() const OMEGA_H_NOEXCEPT {
    return write_.data();
  }
#ifdef OMEGA_H_USE_KOKKOS
  Kokkos::View<const T*> view() const;
#endif
  T get(LO i) const;
  T first() const;
  T last() const;
  OMEGA_H_INLINE bool exists() const OMEGA_H_NOEXCEPT {
    return write_.exists();
  }
  std::string name() const { return write_.name(); }
  OMEGA_H_INLINE T const* begin() const noexcept { return data(); }
  OMEGA_H_INLINE T const* end() const noexcept { return data() + size(); }
};

template <typename T>
Read<T> read(Write<T> a) {
  return Read<T>(a);
}

class Bytes : public Read<Byte> {
 public:
  OMEGA_H_INLINE Bytes() {}
  OMEGA_H_INLINE Bytes(Read<Byte> base) : Read<Byte>(base) {}
  Bytes(Write<Byte> write);
  Bytes(LO size_in, Byte value, std::string const& name = "");
  Bytes(std::initializer_list<Byte> l, std::string const& name = "");
};

class LOs : public Read<LO> {
 public:
  OMEGA_H_INLINE LOs() {}
  OMEGA_H_INLINE LOs(Read<LO> base) : Read<LO>(base) {}
  LOs(Write<LO> write);
  LOs(LO size_in, LO value, std::string const& name = "");
  LOs(LO size_in, LO offset, LO stride, std::string const& name = "");
  LOs(std::initializer_list<LO> l, std::string const& name = "");
};

class GOs : public Read<GO> {
 public:
  OMEGA_H_INLINE GOs() {}
  OMEGA_H_INLINE GOs(Read<GO> base) : Read<GO>(base) {}
  GOs(Write<GO> write);
  GOs(LO size_in, GO value, std::string const& name = "");
  GOs(LO size_in, GO offset, GO stride, std::string const& name = "");
  GOs(std::initializer_list<GO> l, std::string const& name = "");
};

class Reals : public Read<Real> {
 public:
  OMEGA_H_INLINE Reals() {}
  OMEGA_H_INLINE Reals(Read<Real> base) : Read<Real>(base) {}
  Reals(Write<Real> write);
  Reals(LO size_in, Real value, std::string const& name = "");
  Reals(std::initializer_list<Real> l, std::string const& name = "");
};

template <typename T>
class HostRead {
  Read<T> read_;
#if defined(OMEGA_H_USE_KOKKOS)
  typename Kokkos::View<const T*, Kokkos::HostSpace> mirror_;
#elif defined(OMEGA_H_USE_CUDA)
  std::shared_ptr<T[]> mirror_;
#endif
 public:
  typedef T value_type;
  HostRead() = default;
  HostRead(Read<T> read);
  LO size() const;
  inline T const& operator[](LO i) const OMEGA_H_NOEXCEPT {
#ifdef OMEGA_H_CHECK_BOUNDS
    OMEGA_H_CHECK_OP(0, <=, i);
    OMEGA_H_CHECK_OP(i, <, size());
#endif
#ifdef OMEGA_H_USE_KOKKOS
    return mirror_(i);
#else
#ifdef OMEGA_H_USE_CUDA
    return mirror_[i];
#else
    return read_[i];
#endif
#endif
  }
  T const* data() const;
  T get(LO i) const;
  T last() const;
};

template <typename T>
class HostWrite {
  Write<T> write_;
#ifdef OMEGA_H_USE_KOKKOS
  typename Kokkos::View<T*>::HostMirror mirror_;
#elif defined(OMEGA_H_USE_CUDA)
  std::shared_ptr<T[]> mirror_;
#endif
 public:
  typedef T value_type;
  HostWrite() = default;
  HostWrite(LO size_in, std::string const& name = "");
  HostWrite(LO size_in, T offset, T stride, std::string const& name = "");
  HostWrite(Write<T> write_in);
  HostWrite(std::initializer_list<T> l, std::string const& name = "");
  Write<T> write() const;
  LO size() const OMEGA_H_NOEXCEPT;
  inline T& operator[](LO i) const OMEGA_H_NOEXCEPT {
#ifdef OMEGA_H_CHECK_BOUNDS
    OMEGA_H_CHECK_OP(0, <=, i);
    OMEGA_H_CHECK_OP(i, <, size());
#endif
#ifdef OMEGA_H_USE_KOKKOS
    return mirror_(i);
#else
#ifdef OMEGA_H_USE_CUDA
    return mirror_[i];
#else
    return write_[i];
#endif
#endif
  }
  T* data() const;
  OMEGA_H_INLINE bool exists() const OMEGA_H_NOEXCEPT {
    return write_.exists();
  }
  void set(LO i, T value);
  T get(LO i) const;
};

template <typename T>
void fill(Write<T> a, T val);
template <typename T>
void fill_linear(Write<T> a, T offset, T stride);
template <class T>
void copy_into(Read<T> a, Write<T> b);
template <class T>
Write<T> deep_copy(Read<T> a, std::string const& name = "");

/* begin explicit instantiation declarations */
#define OMEGA_H_EXPL_INST_DECL(T)                                              \
  extern template T* nonnull(T*);                                              \
  extern template T const* nonnull(T const*);                                  \
  extern template class Read<T>;                                               \
  extern template class Write<T>;                                              \
  extern template class HostRead<T>;                                           \
  extern template class HostWrite<T>;                                          \
  extern template void fill(Write<T> a, T val);                                \
  extern template void fill_linear(Write<T> a, T, T);                          \
  extern template void copy_into(Read<T> a, Write<T> b);                       \
  extern template Write<T> deep_copy(Read<T> a, std::string const&);
OMEGA_H_EXPL_INST_DECL(I8)
OMEGA_H_EXPL_INST_DECL(I32)
OMEGA_H_EXPL_INST_DECL(I64)
OMEGA_H_EXPL_INST_DECL(Real)
#undef OMEGA_H_EXPL_INST_DECL
/* end explicit instantiation declarations */

}  // end namespace Omega_h

#endif
