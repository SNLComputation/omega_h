#include "Omega_h_array_ops.hpp"

#include "Omega_h_for.hpp"
#include "Omega_h_functors.hpp"
#include "Omega_h_int_iterator.hpp"
#include "Omega_h_reduce.hpp"
#include "Omega_h_dbg.hpp"

namespace Omega_h {

template <typename T>
bool operator==(Read<T> a, Read<T> b) {
  OMEGA_H_CHECK(a.size() == b.size());
  auto const first = IntIterator(0);
  auto const last = IntIterator(a.size());
  auto const init = true;
  auto const op = logical_and<bool>();
  auto transform = OMEGA_H_LAMBDA(LO i)->bool { return a[i] == b[i]; };
  return transform_reduce(first, last, init, op, std::move(transform));
}

template <typename T>
promoted_t<T> get_sum(Read<T> a) {
  using PT = promoted_t<T>;
  return transform_reduce(a.begin(), a.end(), PT(0), plus<PT>(),
      OMEGA_H_LAMBDA(T val)->PT { return PT(val); });
}

template <typename T>
promoted_t<T> get_sum(CommPtr comm, Read<T> a) {
  return comm->allreduce(get_sum(a), OMEGA_H_SUM);
}

template <typename T>
T get_min(Read<T> a) {
  auto const first = IntIterator(0);
  auto const last = IntIterator(a.size());
  auto const init = promoted_t<T>(ArithTraits<T>::max());
  auto const op = minimum<promoted_t<T>>();
  auto transform = OMEGA_H_LAMBDA(LO i)->promoted_t<T> {
    return promoted_t<T>(a[i]);
  };
  auto const r = transform_reduce(first, last, init, op, std::move(transform));
  return T(r);  // see StandinTraits
}

template <typename T>
T get_max(Read<T> a) {
  auto const first = IntIterator(0);
  auto const last = IntIterator(a.size());
  auto const init = promoted_t<T>(ArithTraits<T>::min());
  auto const op = maximum<promoted_t<T>>();
  auto transform = OMEGA_H_LAMBDA(LO i)->promoted_t<T> {
    return promoted_t<T>(a[i]);
  };
  auto const r = transform_reduce(first, last, init, op, std::move(transform));
  return static_cast<T>(r);  // see StandinTraits
}

template <typename T>
T get_min(CommPtr comm, Read<T> a) {
  return comm->allreduce(get_min(a), OMEGA_H_MIN);
}

template <typename T>
T get_max(CommPtr comm, Read<T> a) {
  return comm->allreduce(get_max(a), OMEGA_H_MAX);
}

template <typename T>
MinMax<T> get_minmax(CommPtr comm, Read<T> a) {
  return {get_min(comm, a), get_max(comm, a)};
}

bool are_close(Reals a, Reals b, Real tol, Real floor) {
  OMEGA_H_CHECK(a.size() == b.size());
  auto const first = IntIterator(0);
  auto const last = IntIterator(a.size());
  auto const init = true;
  auto const op = logical_and<bool>();
  auto transform = OMEGA_H_LAMBDA(LO i)->bool {
    return are_close(a[i], b[i], tol, floor);
  };
  auto const res =
      transform_reduce(first, last, init, op, std::move(transform));
  return static_cast<bool>(res);
}

bool are_close_abs(Reals a, Reals b, Real tol) {
  OMEGA_H_CHECK(a.size() == b.size());
  auto const first = IntIterator(0);
  auto const last = IntIterator(a.size());
  auto const init = true;
  auto const op = logical_and<bool>();
  auto transform = OMEGA_H_LAMBDA(LO i)->bool {
    return (std::abs(a[i] - b[i]) <= tol);
  };
  auto const res =
      transform_reduce(first, last, init, op, std::move(transform));
  return static_cast<bool>(res);
}

template <typename T>
Write<T> multiply_each(Read<T> a, Read<T> b, std::string const& name) {
  Write<T> c(a.size(), name);
  if (b.size() == 0) {
    OMEGA_H_CHECK(a.size() == 0);
    return c;
  }
  auto width = divide_no_remainder(a.size(), b.size());
  auto f = OMEGA_H_LAMBDA(LO i) {
    for (Int j = 0; j < width; ++j) {
      c[i * width + j] = a[i * width + j] * b[i];
    }
  };
  parallel_for(b.size(), f, "multiply_each");
  return c;
}

template <typename T>
Read<T> multiply_each_by(Read<T> a, T b) {
  Write<T> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = a[i] * b; };
  parallel_for(a.size(), f, "multiply_each_by");
  return c;
}

template <typename T>
Write<T> divide_each(Read<T> a, Read<T> b, std::string const& name) {
  auto width = divide_no_remainder(a.size(), b.size());
  Write<T> c(a.size(), name);
  auto f = OMEGA_H_LAMBDA(LO i) {
    for (Int j = 0; j < width; ++j) {
      c[i * width + j] = a[i * width + j] / b[i];
    }
  };
  parallel_for(b.size(), f, "divide_each");
  return c;
}

Reals divide_each_maybe_zero(Reals a, Reals b) {
  if (b.size() == 0) {
    OMEGA_H_CHECK(a.size() == 0);
    return a;
  }
  auto width = divide_no_remainder(a.size(), b.size());
  Write<Real> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) {
    if (b[i] != 0.0) {
      for (Int j = 0; j < width; ++j) {
        c[i * width + j] = a[i * width + j] / b[i];
      }
    } else {
      for (Int j = 0; j < width; ++j) {
        OMEGA_H_CHECK(a[i * width + j] == 0.0);
        c[i * width + j] = 0.0;
      }
    }
  };
  parallel_for(b.size(), f, "divide_maybe_zero");
  return c;
}

Reals pow_each(Reals a, Reals b) {
  OMEGA_H_CHECK(a.size() == b.size());
  Write<Real> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = std::pow(a[i], b[i]); };
  parallel_for(a.size(), f, "pow_each");
  return c;
}

template <typename T>
Read<T> divide_each_by(Read<T> a, T b) {
  Write<T> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = a[i] / b; };
  parallel_for(a.size(), f, "divide_each_by");
  return c;
}

template <typename T>
Read<T> add_each(Read<T> a, Read<T> b, std::string const& name) {
  OMEGA_H_CHECK(a.size() == b.size());
  Write<T> c(a.size(), name);
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = a[i] + b[i]; };
  parallel_for(c.size(), f, "add_each");
  return c;
}

template <typename T>
Read<T> subtract_each(Read<T> a, Read<T> b) {
  OMEGA_H_CHECK(a.size() == b.size());
  Write<T> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = a[i] - b[i]; };
  parallel_for(c.size(), f, "subtract_each");
  return c;
}

template <typename T>
Read<T> add_to_each(Read<T> a, T b) {
  Write<T> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = a[i] + b; };
  parallel_for(c.size(), f, "add_to_each");
  return c;
}

template <typename T>
Read<T> subtract_from_each(Read<T> a, T b) {
  Write<T> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = a[i] - b; };
  parallel_for(c.size(), f, "subtract_from_each");
  return c;
}

template <typename T>
Bytes each_geq_to(Read<T> a, T b) {
  Write<I8> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = (a[i] >= b); };
  parallel_for(c.size(), f, "each_geq_to");
  return c;
}

template <typename T>
Bytes each_leq_to(Read<T> a, T b) {
  Write<I8> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = (a[i] <= b); };
  parallel_for(c.size(), f, "each_leq_to");
  return c;
}

template <typename T>
Bytes each_gt(Read<T> a, T b) {
  Write<I8> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = (a[i] > b); };
  parallel_for(c.size(), f, "each_gt");
  return c;
}

template <typename T>
Bytes each_lt(Read<T> a, T b) {
  Write<I8> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = (a[i] < b); };
  parallel_for(c.size(), f, "each_lt");
  return c;
}

template <typename T>
Bytes gt_each(Read<T> a, Read<T> b) {
  OMEGA_H_CHECK(a.size() == b.size());
  Write<I8> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = (a[i] > b[i]); };
  parallel_for(c.size(), f, "gt_each");
  return c;
}

template <typename T>
Bytes lt_each(Read<T> a, Read<T> b) {
  OMEGA_H_CHECK(a.size() == b.size());
  Write<I8> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = (a[i] < b[i]); };
  parallel_for(c.size(), f, "lt_each");
  return c;
}

template <typename T>
Bytes eq_each(Read<T> a, Read<T> b) {
  OMEGA_H_CHECK(a.size() == b.size());
  Write<I8> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = (a[i] == b[i]); };
  parallel_for(c.size(), f, "eq_each");
  return c;
}

template <typename T>
Bytes neq_each(Read<T> a, Read<T> b) {
  OMEGA_H_CHECK(a.size() == b.size());
  Write<I8> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = (a[i] != b[i]); };
  parallel_for("neq_each", c.size(), std::move(f));
  return c;
}

template <typename T>
Bytes geq_each(Read<T> a, Read<T> b) {
  OMEGA_H_CHECK(a.size() == b.size());
  Write<I8> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = (a[i] >= b[i]); };
  parallel_for(c.size(), f, "geq_each");
  return c;
}

template <typename T>
Read<T> min_each(Read<T> a, Read<T> b) {
  OMEGA_H_CHECK(a.size() == b.size());
  Write<T> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = min2(a[i], b[i]); };
  parallel_for(c.size(), f, "min_each");
  return c;
}

template <typename T>
Read<T> max_each(Read<T> a, Read<T> b) {
  OMEGA_H_CHECK(a.size() == b.size());
  Write<T> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = max2(a[i], b[i]); };
  parallel_for(c.size(), f, "max_each");
  return c;
}

template <typename T>
Read<T> ternary_each(Bytes cond, Read<T> a, Read<T> b) {
  OMEGA_H_CHECK(a.size() == b.size());
  auto width = divide_no_remainder(a.size(), cond.size());
  Write<T> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = cond[i / width] ? a[i] : b[i]; };
  parallel_for(c.size(), f, "ternary_each");
  return c;
}

template <typename T>
Read<T> each_max_with(Read<T> a, T b) {
  Write<T> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = max2(a[i], b); };
  parallel_for(c.size(), f, "each_max_with");
  return c;
}

template <typename T>
Bytes each_neq_to(Read<T> a, T b) {
  Write<I8> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = (a[i] != b); };
  parallel_for(c.size(), f, "each_neq_to");
  return c;
}

template <typename T>
Bytes each_eq(Read<T> a, Read<T> b) {
  OMEGA_H_CHECK(a.size() == b.size());
  Write<I8> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = (a[i] == b[i]); };
  parallel_for(c.size(), f, "each_eq");
  return c;
}

template <typename T>
Bytes each_eq_to(Read<T> a, T b) {
  Write<I8> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = (a[i] == b); };
  parallel_for(c.size(), f, "each_eq_to");
  return c;
}

Bytes land_each(Bytes a, Bytes b) {
  OMEGA_H_CHECK(a.size() == b.size());
  Write<I8> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = (a[i] && b[i]); };
  parallel_for(c.size(), f, "land_each");
  return c;
}

Bytes lor_each(Bytes a, Bytes b) {
  OMEGA_H_CHECK(a.size() == b.size());
  Write<I8> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = (a[i] || b[i]); };
  parallel_for(c.size(), f, "lor_each");
  return c;
}

Bytes bit_or_each(Bytes a, Bytes b) {
  OMEGA_H_CHECK(a.size() == b.size());
  Write<I8> c(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { c[i] = (a[i] | b[i]); };
  parallel_for(c.size(), f, "bit_or_each");
  return c;
}

Bytes bit_neg_each(Bytes a) {
  Write<I8> b(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { b[i] = ~(a[i]); };
  parallel_for(a.size(), f, "bit_neg_each");
  return b;
}

Read<Real> fabs_each(Read<Real> a) {
  Write<Real> b(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { b[i] = std::abs(a[i]); };
  parallel_for(a.size(), f, "fabs_each");
  return b;
}

template <typename T>
Read<T> get_component(Read<T> a, Int ncomps, Int comp) {
  Write<T> b(divide_no_remainder(a.size(), ncomps));
  auto f = OMEGA_H_LAMBDA(LO i) { b[i] = a[i * ncomps + comp]; };
  parallel_for(b.size(), f, "get_component");
  return b;
}

template <typename T>
void set_component(Write<T> out, Read<T> a, Int ncomps, Int comp) {
  auto f = OMEGA_H_LAMBDA(LO i) { out[i * ncomps + comp] = a[i]; };
  parallel_for(a.size(), f, "set_component");
}

template <typename T>
LO find_last(Read<T> array, T value) {
  auto const first = IntIterator(0);
  auto const last = IntIterator(array.size());
  auto const init = -1;
  auto const op = maximum<LO>();
  auto transform = OMEGA_H_LAMBDA(LO i)->LO {
    if (array[i] == value)
      return i;
    else
      return -1;
  };
  return transform_reduce(first, last, init, op, std::move(transform));
}

template <typename T>
bool is_sorted(Read<T> a) {
  if (a.size() < 2) return true;
  auto const first = IntIterator(0);
  auto const last = IntIterator(a.size() - 1);
  auto const init = true;
  auto const op = logical_and<bool>();
  auto transform = OMEGA_H_LAMBDA(LO i)->bool { return a[i] <= a[i + 1]; };
  return transform_reduce(first, last, init, op, std::move(transform));
}

template <typename T>
Read<T> interleave(std::vector<Read<T>> arrays) {
  if (arrays.empty()) return Read<T>();
  auto narrays = LO(arrays.size());
  auto array_size = arrays.front().size();
  for (auto& array : arrays) OMEGA_H_CHECK(array.size() == array_size);
  auto out_size = narrays * array_size;
  auto out = Write<T>(out_size);
  for (LO i = 0; i < narrays; ++i) {
    auto array = arrays[std::size_t(i)];
    auto f = OMEGA_H_LAMBDA(LO j) { out[j * narrays + i] = array[j]; };
    parallel_for(array_size, f, "interleave");
  }
  return out;
}

template <typename T>
Read<T> coalesce(std::vector<Read<T>> arrays) {
  if (arrays.empty()) return Read<T>();
  std::vector<LO> offsets(arrays.size() + 1);
  OMEGA_H_CHECK(offsets.data() != nullptr);
  offsets[0] = 0;
  for (std::size_t i = 1; i <= arrays.size(); ++i) {
    offsets[i] = offsets[i - 1] + arrays[i].size();
  }
  auto out_size = offsets[arrays.size()];
  auto out = Write<T>(out_size);
  for (std::size_t i = 0; i < arrays.size(); ++i) {
    auto array = arrays[std::size_t(i)];
    auto offset = offsets[i];
    auto f = OMEGA_H_LAMBDA(LO j) { out[offset + j] = array[j]; };
    parallel_for(array.size(), f, "coalesce");
  }
  return out;
}

/* A reproducible sum of floating-point values.
   this operation is one of the key places where
   a program's output begins to depend on parallel
   partitioning and traversal order, because
   floating-point values do not produce the same
   sum when added in a different order.

   IEEE 754 64-bit floating point format is assumed,
   which has 52 bits in the fraction.

   The idea here is to add the numbers as fixed-point values.
   max_exponent() finds the largest exponent (e) such that
   all values are (<= 2^(e)).
   We then use the value (2^(e - 52)) as the unit, and sum all
   values as integers in that unit.
   This is guaranteed to be at least as accurate as the
   worst-case ordering of the values, i.e. being added
   in order of decreasing magnitude.

   If we used a 64-bit integer type, we would only be
   able to reliably add up to (2^12 = 4096) values
   (64 - 52 = 12).
   Thus we use a 128-bit integer type.
   This allows us to reliably add up to (2^76 > 10^22)
   values. By comparison, supercomputers today
   support a maximum of one million MPI ranks (10^6)
   and each rank typically can't hold more than
   one billion values (10^9), for a total of (10^15) values.
*/

int max_exponent(Reals a) {
  if (a.size() == 0) {
    return 0;
  }
  auto const first = IntIterator(0);
  auto const last = IntIterator(a.size());
  auto const init = ArithTraits<int>::min();
  auto const op = maximum<int>();
  auto transform = OMEGA_H_LAMBDA(LO i)->int {
    int expo;
    std::frexp(a[i], &expo);
    return expo;
  };
  return transform_reduce(first, last, init, op, std::move(transform));
}

struct Int128Plus {
  OMEGA_H_INLINE Int128 operator()(Int128 a, Int128 b) const { return a + b; }
};

Int128 int128_sum(Reals const a, double const unit) {
  if (a.size() == 0) {
    return Int128(0);
  }
  auto const first = IntIterator(0);
  auto const last = IntIterator(a.size());
  auto const init = Int128(0);
  auto const op = Int128Plus();
  auto transform = OMEGA_H_LAMBDA(LO i)->Int128 {
    return Int128::from_double(a[i], unit);
  };
  return transform_reduce(first, last, init, op, std::move(transform));
}

Real repro_sum(Reals a) {
  if (a.size() == 0) {
    return 0.0;
  }
  begin_code("repro_sum");
  int expo = max_exponent(a);
  double unit = exp2(double(expo - MANTISSA_BITS));
  Int128 fixpt_sum = int128_sum(a, unit);
  end_code();
  return fixpt_sum.to_double(unit);
}

Real repro_sum(CommPtr comm, Reals a) {
  begin_code("repro_sum(comm)");
  int expo = comm->allreduce(max_exponent(a), OMEGA_H_MAX);
  double unit = exp2(double(expo - MANTISSA_BITS));
  Int128 fixpt_sum = int128_sum(a, unit);
  fixpt_sum = comm->add_int128(fixpt_sum);
  end_code();
  return fixpt_sum.to_double(unit);
}

void repro_sum(CommPtr comm, Reals a, Int ncomps, Real result[]) {
  for (Int comp = 0; comp < ncomps; ++comp) {
    result[comp] = repro_sum(comm, get_component(a, ncomps, comp));
  }
}

Reals interpolate_between(Reals a, Reals b, Real t) {
  OMEGA_H_CHECK(a.size() == b.size());
  auto n = a.size();
  auto out = Write<Real>(n);
  auto f = OMEGA_H_LAMBDA(LO i) { out[i] = a[i] * (1.0 - t) + b[i] * t; };
  parallel_for(n, f, "interpolate_between");
  return out;
}

Reals invert_each(Reals a) {
  auto out = Write<Real>(a.size());
  auto f = OMEGA_H_LAMBDA(LO i) { out[i] = 1.0 / a[i]; };
  parallel_for(a.size(), f, "invert_each");
  return out;
}

template <typename Tout, typename Tin>
Read<Tout> array_cast(Read<Tin> in) {
  auto out = Write<Tout>(in.size());
  auto f = OMEGA_H_LAMBDA(LO i) { out[i] = static_cast<Tout>(in[i]); };
  parallel_for(in.size(), f, "array_cast");
  return out;
}

#define INST(T)                                                                \
  template bool operator==(Read<T> a, Read<T> b);                              \
  template promoted_t<T> get_sum(Read<T> a);                                   \
  template T get_min(Read<T> a);                                               \
  template T get_max(Read<T> a);                                               \
  template promoted_t<T> get_sum(CommPtr comm, Read<T> a);                     \
  template T get_min(CommPtr comm, Read<T> a);                                 \
  template T get_max(CommPtr comm, Read<T> a);                                 \
  template MinMax<T> get_minmax(CommPtr comm, Read<T> a);                      \
  template Write<T> multiply_each(Read<T> a, Read<T> b, std::string const&);   \
  template Read<T> multiply_each_by(Read<T> a, T b);                           \
  template Read<T> divide_each_by(Read<T> x, T b);                             \
  template Write<T> divide_each(Read<T> a, Read<T> b, std::string const&);     \
  template Read<T> add_each(Read<T> a, Read<T> b, std::string const&);         \
  template Read<T> subtract_each(Read<T> a, Read<T> b);                        \
  template Read<T> min_each(Read<T> a, Read<T> b);                             \
  template Read<T> max_each(Read<T> a, Read<T> b);                             \
  template Read<T> ternary_each(Bytes cond, Read<T> a, Read<T> b);             \
  template Read<T> each_max_with(Read<T> a, T b);                              \
  template Read<T> add_to_each(Read<T> a, T b);                                \
  template Read<T> subtract_from_each(Read<T> a, T b);                         \
  template Bytes each_geq_to(Read<T> a, T b);                                  \
  template Bytes each_leq_to(Read<T> a, T b);                                  \
  template Bytes each_gt(Read<T> a, T b);                                      \
  template Bytes each_lt(Read<T> a, T b);                                      \
  template Bytes each_neq_to(Read<T> a, T b);                                  \
  template Bytes each_eq(Read<T> a, Read<T> b);                                \
  template Bytes each_eq_to(Read<T> a, T b);                                   \
  template Bytes gt_each(Read<T> a, Read<T> b);                                \
  template Bytes lt_each(Read<T> a, Read<T> b);                                \
  template Bytes eq_each(Read<T> a, Read<T> b);                                \
  template Bytes neq_each(Read<T> a, Read<T> b);                               \
  template Read<T> get_component(Read<T> a, Int ncomps, Int comp);             \
  template void set_component(Write<T> out, Read<T> a, Int ncomps, Int comp);  \
  template LO find_last(Read<T> array, T value);                               \
  template bool is_sorted(Read<T> a);                                          \
  template Read<T> interleave(std::vector<Read<T>> arrays);                    \
  template Read<T> coalesce(std::vector<Read<T>> arrays);

INST(I8)
INST(I32)
INST(I64)
INST(Real)
#undef INST

template Read<Real> array_cast(Read<I32>);
template Read<I32> array_cast(Read<I8>);

}  // end namespace Omega_h
