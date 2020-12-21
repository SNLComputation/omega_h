#ifndef OMEGA_H_ARRAY_OPS_HPP
#define OMEGA_H_ARRAY_OPS_HPP

#include <vector>

#include <Omega_h_comm.hpp>
#include <Omega_h_scalar.hpp>

namespace Omega_h {

template <class T>
bool operator==(Read<T> a, Read<T> b);

template <typename T>
promoted_t<T> get_sum(Read<T> a);
template <typename T>
T get_min(Read<T> a);
template <typename T>
T get_max(Read<T> a);

template <typename T>
promoted_t<T> get_sum(CommPtr comm, Read<T> a);
template <typename T>
T get_min(CommPtr comm, Read<T> a);
template <typename T>
T get_max(CommPtr comm, Read<T> a);

template <typename T>
struct MinMax {
  T min;
  T max;
};

template <typename T>
MinMax<T> get_minmax(CommPtr comm, Read<T> a);

bool are_close(Reals a, Reals b, Real tol = EPSILON, Real floor = EPSILON);
bool are_close_abs(Reals a, Reals b, Real tol = EPSILON);

template <typename T>
Bytes each_eq(Read<T> a, Read<T> b);
template <typename T>
Bytes each_eq_to(Read<T> a, T b);

/* "a" may be larger than "b" by some integer factor */
template <typename T>
Write<T> multiply_each(Read<T> a, Read<T> b, std::string const& name = "");
template <typename T>
Read<T> multiply_each_by(Read<T> a, T b);
template <typename T>
Write<T> divide_each(Read<T> a, Read<T> b, std::string const& name = "");
Reals divide_each_maybe_zero(Reals a, Reals b);
template <typename T>
Write<T> concat(Read<T> a, Read<T> b, std::string const& name = "");
Reals pow_each(Reals a, Reals b);
template <typename T>
Read<T> divide_each_by(Read<T> a, T b);
template <typename T>
Read<T> add_each(Read<T> a, Read<T> b, std::string const& name = "");
template <typename T>
Read<T> subtract_each(Read<T> a, Read<T> b);
template <typename T>
Read<T> add_to_each(Read<T> a, T b);
template <typename T>
Read<T> subtract_from_each(Read<T> a, T b);
template <typename T>
Bytes each_geq_to(Read<T> a, T b);
template <typename T>
Bytes each_leq_to(Read<T> a, T b);
template <typename T>
Bytes each_gt(Read<T> a, T b);
template <typename T>
Bytes each_lt(Read<T> a, T b);
template <typename T>
Bytes gt_each(Read<T> a, Read<T> b);
template <typename T>
Bytes lt_each(Read<T> a, Read<T> b);
template <typename T>
Bytes eq_each(Read<T> a, Read<T> b);
template <typename T>
Bytes neq_each(Read<T> a, Read<T> b);
template <typename T>
Bytes each_neq_to(Read<T> a, T b);
template <typename T>
Read<T> min_each(Read<T> a, Read<T> b);
template <typename T>
Read<T> max_each(Read<T> a, Read<T> b);
template <typename T>
Read<T> each_max_with(Read<T> a, T b);
Bytes land_each(Bytes a, Bytes b);
Bytes lor_each(Bytes a, Bytes b);
Bytes bit_or_each(Bytes a, Bytes b);
Bytes bit_neg_each(Bytes a);
Read<Real> fabs_each(Read<Real> a);

template <typename T>
Read<T> ternary_each(Bytes cond, Read<T> a, Read<T> b);

template <typename T>
Read<T> get_component(Read<T> a, Int ncomps, Int comp);

template <typename T>
void set_component(Write<T> out, Read<T> a, Int ncomps, Int comp);

template <typename T>
LO find_last(Read<T> array, T value);

template <typename T>
bool is_sorted(Read<T> array);

template <typename T>
Read<T> interleave(std::vector<Read<T>> arrays);

template <typename T>
Read<T> coalesce(std::vector<Read<T>> arrays);

Real repro_sum(Reals a);
Real repro_sum(CommPtr comm, Reals a);
void repro_sum(CommPtr comm, Reals a, Int ncomps, Real result[]);

Reals interpolate_between(Reals a, Reals b, Real t);
Reals invert_each(Reals a);

template <typename Tout, typename Tin>
Read<Tout> array_cast(Read<Tin> in);

#define OMEGA_H_EXPL_INST_DECL(T)                                              \
  extern template bool operator==(Read<T> a, Read<T> b);                       \
  extern template promoted_t<T> get_sum(Read<T> a);                            \
  extern template T get_min(Read<T> a);                                        \
  extern template T get_max(Read<T> a);                                        \
  extern template promoted_t<T> get_sum(CommPtr comm, Read<T> a);              \
  extern template T get_min(CommPtr comm, Read<T> a);                          \
  extern template T get_max(CommPtr comm, Read<T> a);                          \
  extern template MinMax<T> get_minmax(CommPtr comm, Read<T> a);               \
  extern template Write<T> multiply_each(                                      \
      Read<T> a, Read<T> b, std::string const&);                               \
  extern template Write<T> divide_each(                                        \
      Read<T> a, Read<T> b, std::string const&);                               \
  extern template Write<T> concat(                                        \
      Read<T> a, Read<T> b, std::string const&);                               \
  extern template Read<T> add_each(Read<T> a, Read<T> b, std::string const&);  \
  extern template Read<T> subtract_each(Read<T> a, Read<T> b);                 \
  extern template Read<T> add_to_each(Read<T> a, T b);                         \
  extern template Read<T> subtract_from_each(Read<T> a, T b);                  \
  extern template Bytes each_geq_to(Read<T> a, T b);                           \
  extern template Bytes each_leq_to(Read<T> a, T b);                           \
  extern template Bytes each_gt(Read<T> a, T b);                               \
  extern template Bytes each_lt(Read<T> a, T b);                               \
  extern template Bytes each_neq_to(Read<T> a, T b);                           \
  extern template Bytes each_eq(Read<T> a, Read<T> b);                         \
  extern template Bytes each_eq_to(Read<T> a, T b);                            \
  extern template Read<T> multiply_each_by(Read<T> a, T b);                    \
  extern template Read<T> divide_each_by(Read<T> a, T b);                      \
  extern template Read<T> min_each(Read<T> a, Read<T> b);                      \
  extern template Read<T> max_each(Read<T> a, Read<T> b);                      \
  extern template Read<T> ternary_each(Bytes cond, Read<T> a, Read<T> b);      \
  extern template Read<T> each_max_with(Read<T> a, T b);                       \
  extern template Bytes gt_each(Read<T> a, Read<T> b);                         \
  extern template Bytes lt_each(Read<T> a, Read<T> b);                         \
  extern template Bytes eq_each(Read<T> a, Read<T> b);                         \
  extern template Bytes neq_each(Read<T> a, Read<T> b);                        \
  extern template Read<T> get_component(Read<T> a, Int ncomps, Int comp);      \
  extern template void set_component(                                          \
      Write<T> out, Read<T> a, Int ncomps, Int comp);                          \
  extern template LO find_last(Read<T> array, T value);                        \
  extern template bool is_sorted(Read<T> array);                               \
  extern template Read<T> interleave(std::vector<Read<T>> arrays);             \
  extern template Read<T> coalesce(std::vector<Read<T>> arrays);
OMEGA_H_EXPL_INST_DECL(I8)
OMEGA_H_EXPL_INST_DECL(I32)
OMEGA_H_EXPL_INST_DECL(I64)
OMEGA_H_EXPL_INST_DECL(Real)
#undef OMEGA_H_EXPL_INST_DECL

extern template Read<Real> array_cast(Read<I32>);
extern template Read<I32> array_cast(Read<I8>);

}  // end namespace Omega_h

#endif
