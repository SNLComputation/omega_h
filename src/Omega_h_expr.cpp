#include <sstream>

#include <Omega_h_array_ops.hpp>
#include <Omega_h_expr.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_math_lang.hpp>
#include <Omega_h_matrix.hpp>
#include <Omega_h_vector.hpp>
#include <Omega_h_dbg.hpp>

#include <sstream>

namespace Omega_h {

namespace {

void promote_bool(LO size, any& x) { x = Bytes(size, Byte(any_cast<bool>(x))); }

void promote_bools(LO size, any& lhs, any& rhs) {
  if (lhs.type() == typeid(bool) && rhs.type() == typeid(Bytes)) {
    OMEGA_H_CHECK(any_cast<Bytes>(rhs).size() == size);
    promote_bool(size, lhs);
  }
  if (lhs.type() == typeid(Bytes) && rhs.type() == typeid(bool)) {
    OMEGA_H_CHECK(any_cast<Bytes>(lhs).size() == size);
    promote_bool(size, rhs);
  }
}

void promote_scalar(LO size, any& x) { x = Reals(size, any_cast<Real>(x)); }

void promote_scalars(LO size, any& lhs, any& rhs) {
  if (lhs.type() == typeid(Real) && rhs.type() == typeid(Reals)) {
    promote_scalar(size, lhs);
  }
  if (lhs.type() == typeid(Reals) && rhs.type() == typeid(Real)) {
    promote_scalar(size, rhs);
  }
}

template <Int dim>
void promote_vector(LO size, any& x) {
  x = repeat_vector(size, any_cast<Vector<dim>>(x));
}

template <Int dim>
void promote_vectors(LO size, any& lhs, any& rhs) {
  if (lhs.type() == typeid(Vector<dim>) && rhs.type() == typeid(Reals)) {
    promote_vector<dim>(size, lhs);
  }
  if (lhs.type() == typeid(Reals) && rhs.type() == typeid(Vector<dim>)) {
    promote_vector<dim>(size, rhs);
  }
}

template <Int dim>
void promote_matrix(LO size, any& x) {
  x = repeat_matrix(size, any_cast<Tensor<dim>>(x));
}

template <Int dim>
void promote_symm(LO size, any& x) {
  x = repeat_symm(size, vector2symm(any_cast<Vector<symm_ncomps(dim)>>(x)));
}

template <Int dim>
void promote_matrices(LO size, any& lhs, any& rhs) {
  if (lhs.type() == typeid(Tensor<dim>) && rhs.type() == typeid(Reals)) {
    promote_matrix<dim>(size, lhs);
  }
  if (lhs.type() == typeid(Reals) && rhs.type() == typeid(Tensor<dim>)) {
    promote_matrix<dim>(size, rhs);
  }
}

template <Int dim>
void promote(LO size, any& lhs, any& rhs) {
  if (lhs.type() == rhs.type()) return;
  promote_bools(size, lhs, rhs);
  promote_scalars(size, lhs, rhs);
  promote_vectors<dim>(size, lhs, rhs);
  promote_matrices<dim>(size, lhs, rhs);
}

void promote(LO size, Int dim, any& lhs, any& rhs) {
  ScopedTimer("promote(size, dim, lhs, rhs)");
  if (dim == 3)
    promote<3>(size, lhs, rhs);
  else if (dim == 2)
    promote<2>(size, lhs, rhs);
  else if (dim == 1)
    promote<1>(size, lhs, rhs);
  else {
    OMEGA_H_NORETURN();
  }
}

template <Int dim>
void promote(LO size, any& x) {
  if (x.type() == typeid(bool)) {
    promote_bool(size, x);
  } else if (x.type() == typeid(Real)) {
    promote_scalar(size, x);
  } else if (x.type() == typeid(Vector<dim>)) {
    promote_vector<dim>(size, x);
  } else if (x.type() == typeid(Tensor<dim>)) {
    promote_matrix<dim>(size, x);
  } else if (x.type() == typeid(Vector<symm_ncomps(dim)>)) {
    promote_symm<dim>(size, x);
  } else if (x.type() == typeid(Reals)) {
    return;
  } else {
    std::string msg;
    msg += "Unexpected type ";
    msg += x.type().name();
    msg += " being promoted";
    throw ParserFail(msg);
  }
}

void promote(LO size, Int dim, any& x) {
  ScopedTimer("promote(size, dim, any)");
  if (dim == 3)
    promote<3>(size, x);
  else if (dim == 2)
    promote<2>(size, x);
  else if (dim == 1)
    promote<1>(size, x);
  else
    OMEGA_H_NORETURN();
}

any ternary(LO size, Int dim, any const& cond, any& true_val, any& false_val) {
  OMEGA_H_TIME_FUNCTION;
  if (cond.type() == typeid(bool)) {
    if (true_val.type() == typeid(Real)) {
      return any_cast<bool>(cond) ? any_cast<Real>(true_val)
                                  : any_cast<Real>(false_val);
    } else if (true_val.type() == typeid(Reals)) {
      return any_cast<bool>(cond) ? any_cast<Reals>(true_val)
                                  : any_cast<Reals>(false_val);
    } else {
      throw ParserFail("Invalid true value type in ternary operator");
    }
  } else if (cond.type() == typeid(Bytes)) {
    promote(size, dim, true_val);
    promote(size, dim, false_val);
    return Reals(ternary_each(any_cast<Bytes>(cond), any_cast<Reals>(true_val),
        any_cast<Reals>(false_val)));
  } else {
    throw ParserFail("Invalid condition value type in ternary operator");
  }
}

any eval_or(any& lhs, any& rhs) {
  if (lhs.type() == typeid(bool)) {
    return any_cast<bool>(lhs) || any_cast<bool>(rhs);
  } else if (lhs.type() == typeid(Bytes)) {
    return lor_each(any_cast<Bytes>(lhs), any_cast<Bytes>(rhs));
  } else {
    throw ParserFail("Invalid operand types to || operator");
  }
}

any eval_and(any& lhs, any& rhs) {
  if (lhs.type() == typeid(bool)) {
    return any_cast<bool>(lhs) && any_cast<bool>(rhs);
  } else if (lhs.type() == typeid(Bytes)) {
    return land_each(any_cast<Bytes>(lhs), any_cast<Bytes>(rhs));
  } else {
    throw ParserFail("Invalid operand types to && operator");
  }
}

any gt(any& lhs, any& rhs) {
  if (lhs.type() == typeid(Real)) {
    return any_cast<Real>(lhs) > any_cast<Real>(rhs);
  } else if (lhs.type() == typeid(Reals)) {
    return gt_each(any_cast<Reals>(lhs), any_cast<Reals>(rhs));
  } else {
    std::stringstream ss;
    ss << "Invalid operand types to > operator: " << lhs.type().name() << ", "
       << rhs.type().name() << '\n';
    throw ParserFail(ss.str());
  }
}

any lt(any& lhs, any& rhs) {
  if (lhs.type() == typeid(Real)) {
    return any_cast<Real>(lhs) < any_cast<Real>(rhs);
  } else if (lhs.type() == typeid(Reals)) {
    return lt_each(any_cast<Reals>(lhs), any_cast<Reals>(rhs));
  } else {
    throw ParserFail("Invalid operand types to < operator");
  }
}

any eq(any& lhs, any& rhs) {
  if (lhs.type() == typeid(Real)) {
    return any_cast<Real>(lhs) == any_cast<Real>(rhs);
  } else if (lhs.type() == typeid(Reals)) {
    return eq_each(any_cast<Reals>(lhs), any_cast<Reals>(rhs));
  } else {
    throw ParserFail("Invalid operand types to == operator");
  }
}

template <Int dim>
any add(any& lhs, any& rhs) {
  if (lhs.type() == typeid(Real)) {
    return any_cast<Real>(lhs) + any_cast<Real>(rhs);
  } else if (lhs.type() == typeid(Vector<dim>)) {
    return any_cast<Vector<dim>>(lhs) + any_cast<Vector<dim>>(rhs);
  } else if (lhs.type() == typeid(Tensor<dim>)) {
    return any_cast<Tensor<dim>>(lhs) + any_cast<Tensor<dim>>(rhs);
  } else if (lhs.type() == typeid(Reals)) {
    return Reals(add_each(any_cast<Reals>(lhs), any_cast<Reals>(rhs)));
  } else {
    throw ParserFail("Invalid operand types to + operator");
  }
}

any add(Int dim, any& lhs, any& rhs) {
  if (dim == 3) return add<3>(lhs, rhs);
  if (dim == 2)
    return add<2>(lhs, rhs);
  else
    return add<1>(lhs, rhs);
}

template <Int dim>
any sub(any& lhs, any& rhs) {
  if (lhs.type() == typeid(Real)) {
    return any_cast<Real>(lhs) - any_cast<Real>(rhs);
  } else if (lhs.type() == typeid(Vector<dim>)) {
    return any_cast<Vector<dim>>(lhs) - any_cast<Vector<dim>>(rhs);
  } else if (lhs.type() == typeid(Tensor<dim>)) {
    return any_cast<Tensor<dim>>(lhs) - any_cast<Tensor<dim>>(rhs);
  } else if (lhs.type() == typeid(Reals)) {
    return Reals(subtract_each(any_cast<Reals>(lhs), any_cast<Reals>(rhs)));
  } else {
    throw ParserFail("Invalid operand types to + operator");
  }
}

any sub(Int dim, any& lhs, any& rhs) {
  if (dim == 3) return sub<3>(lhs, rhs);
  if (dim == 2)
    return sub<2>(lhs, rhs);
  else
    return sub<1>(lhs, rhs);
}

template <Int dim>
any mul(LO size, any& lhs, any& rhs) {
  if (lhs.type() == typeid(Real) && rhs.type() == typeid(Real)) {
    return any_cast<Real>(lhs) * any_cast<Real>(rhs);
    /* begin multiply non-scalar by scalar (commutative) */
  } else if (lhs.type() == typeid(Vector<dim>) && rhs.type() == typeid(Real)) {
    return any_cast<Vector<dim>>(lhs) * any_cast<Real>(rhs);
  } else if (lhs.type() == typeid(Real) && rhs.type() == typeid(Vector<dim>)) {
    return any_cast<Real>(lhs) * any_cast<Vector<dim>>(rhs);
  } else if (lhs.type() == typeid(Tensor<dim>) && rhs.type() == typeid(Real)) {
    return any_cast<Tensor<dim>>(lhs) * any_cast<Real>(rhs);
  } else if (lhs.type() == typeid(Real) && rhs.type() == typeid(Tensor<dim>)) {
    return any_cast<Real>(lhs) * any_cast<Tensor<dim>>(rhs);
    /* dot product */
  } else if (lhs.type() == typeid(Vector<dim>) &&
             rhs.type() == typeid(Vector<dim>)) {
    return any_cast<Vector<dim>>(lhs) * any_cast<Vector<dim>>(rhs);
    /* matrix * vector (non-commutative) */
  } else if (lhs.type() == typeid(Tensor<dim>) &&
             rhs.type() == typeid(Vector<dim>)) {
    return any_cast<Tensor<dim>>(lhs) * any_cast<Vector<dim>>(rhs);
    /* matrix * matrix (non-commutative) */
  } else if (lhs.type() == typeid(Tensor<dim>) &&
             rhs.type() == typeid(Tensor<dim>)) {
    return any_cast<Tensor<dim>>(lhs) * any_cast<Tensor<dim>>(rhs);
  } else if (lhs.type() == typeid(Reals) && rhs.type() == typeid(Reals)) {
    auto& lhs_vals = any_cast<Reals&>(lhs);
    auto& rhs_vals = any_cast<Reals&>(rhs);
    if (rhs_vals.size() == size) {  // RHS is scalars
      return Reals(multiply_each(lhs_vals, rhs_vals));
    } else if (lhs_vals.size() == size) {  // LHS is scalars
      return Reals(multiply_each(rhs_vals, lhs_vals));
    } else if (lhs_vals.size() == size * dim &&
               rhs_vals.size() == size * dim) {  // dot products
      return dot_vectors(lhs_vals, rhs_vals, dim);
    } else if (lhs_vals.size() == size * matrix_ncomps(dim, dim) &&
               rhs_vals.size() == size * dim) {  // matrices * vectors
      return matrices_times_vectors(lhs_vals, rhs_vals, dim);
    } else if (lhs_vals.size() == size * matrix_ncomps(dim, dim) &&
               rhs_vals.size() == size * matrix_ncomps(dim, dim)) {
      return matrices_times_matrices(lhs_vals, rhs_vals, dim);
    } else {
      throw ParserFail("Unexpected array size in * operator");
    }
  } else {
    throw ParserFail("Invalid operand types to * operator");
  }
}

any mul(LO size, Int dim, any& lhs, any& rhs) {
  if (dim == 3) return mul<3>(size, lhs, rhs);
  if (dim == 2)
    return mul<2>(size, lhs, rhs);
  else
    return mul<1>(size, lhs, rhs);
}

template <Int dim>
any div(any& lhs, any& rhs) {
  if (rhs.type() == typeid(Reals)) {
    return Reals(
        divide_each_maybe_zero(any_cast<Reals>(lhs), any_cast<Reals>(rhs)));
  } else if (rhs.type() == typeid(Real)) {
    if (lhs.type() == typeid(Real)) {
      return any_cast<Real>(lhs) / any_cast<Real>(rhs);
    } else if (lhs.type() == typeid(Vector<dim>)) {
      return any_cast<Vector<dim>>(lhs) / any_cast<Real>(rhs);
    } else if (lhs.type() == typeid(Tensor<dim>)) {
      return any_cast<Tensor<dim>>(lhs) / any_cast<Real>(rhs);
    } else {
      throw ParserFail("Invalid left operand type in / operator");
    }
  } else {
    throw ParserFail("Invalid right operand type in / operator");
  }
}

any div(Int dim, any& lhs, any& rhs) {
  if (dim == 3) return div<3>(lhs, rhs);
  if (dim == 2) return div<2>(lhs, rhs);
  return div<1>(lhs, rhs);
}

template <Int dim>
any eval_pow(any& lhs, any& rhs) {
  if (rhs.type() == typeid(Reals)) {
    return Reals(pow_each(any_cast<Reals>(lhs), any_cast<Reals>(rhs)));
  } else if (rhs.type() == typeid(Real)) {
    return std::pow(any_cast<Real>(lhs), any_cast<Real>(rhs));
  } else {
    throw ParserFail("Invalid right operand type in ^ operator");
  }
}

any eval_pow(Int dim, any& lhs, any& rhs) {
  if (dim == 3) return eval_pow<3>(lhs, rhs);
  if (dim == 2) return eval_pow<2>(lhs, rhs);
  return eval_pow<1>(lhs, rhs);
}

template <Int dim>
any neg_dim(any const& val) {
  if (val.type() == typeid(Real)) {
    return -any_cast<Real>(val);
  } else if (val.type() == typeid(Vector<dim>)) {
    return -any_cast<Vector<dim>>(val);
  } else if (val.type() == typeid(Tensor<dim>)) {
    return -any_cast<Tensor<dim>>(val);
  } else if (val.type() == typeid(Reals)) {
    return Reals(multiply_each_by(any_cast<Reals>(val), -1.0));
  } else {
    throw ParserFail("Invalid operand type to negation operator");
  }
}

any neg(Int dim, any const& val) {
  if (dim == 3) return neg_dim<3>(val);
  if (dim == 2) return neg_dim<2>(val);
  return neg_dim<1>(val);
}

template <Int dim>
any access(LO size, any& var, ExprReader::Args& args) {
  auto i = static_cast<Int>(any_cast<Real>(args.at(0)));
  auto j =
      args.size() > 1 ? static_cast<Int>(any_cast<Real>(args.at(0))) : Int(-1);
  if (var.type() == typeid(Vector<dim>)) {
    return (any_cast<Vector<dim>>(var))(i);
  } else if (var.type() == typeid(Tensor<dim>)) {
    return (any_cast<Tensor<dim>>(var))(i, j);
  } else if (var.type() == typeid(Reals)) {
    auto array = any_cast<Reals>(var);
    if (array.size() == size * dim) {
      return Reals(get_component(array, dim, i));
    } else if (array.size() == size * matrix_ncomps(dim, dim)) {
      return Reals(get_component(array, matrix_ncomps(dim, dim), j * dim + i));
    } else {
      std::stringstream ss;
      ss << "Unexpected array size " << array.size() << " in access operator\n";
      ss << "Value count is " << size << " dimension is " << dim << '\n';
      throw ParserFail(ss.str());
    }
  } else {
    throw ParserFail("Unexpected variable type in access operator\n");
  }
}

any access(LO size, Int dim, any& var, ExprReader::Args& args) {
  if (dim == 3) return access<3>(size, var, args);
  if (dim == 2) return access<2>(size, var, args);
  return access<1>(size, var, args);
}

template <Int dim>
any make_vector(ExprReader::Args& args) {
  auto v = zero_vector<dim>();
  Int i;
  for (i = 0; i < Int(args.size()); ++i) {
    auto& arg = args[std::size_t(i)];
    v[i] = any_cast<Real>(arg);
  }
  for (; i < dim; ++i) {
    v[i] = v[Int(args.size() - 1)];
  }
  return v;
}

any make_vector(LO size, Int dim, ExprReader::Args& args) {
  if (args.size() > std::size_t(dim)) {
    throw ParserFail("Too many arguments to vector()\n");
  }
  bool has_arrays = false;
  for (auto& arg : args)
    if (arg.type() == typeid(Reals)) has_arrays = true;
  if (has_arrays) {
    std::vector<Read<Real>> arrays;
    Int i;
    for (i = 0; i < Int(args.size()); ++i) {
      auto& arg = args[std::size_t(i)];
      promote(size, dim, arg);
      arrays.push_back(any_cast<Reals>(arg));
    }
    for (; i < dim; ++i) {
      arrays.push_back(arrays.back());
    }
    OMEGA_H_CHECK(Int(arrays.size()) == dim);
    return Reals(interleave(arrays));
  } else {
    if (dim == 3) return make_vector<3>(args);
    if (dim == 2) return make_vector<2>(args);
    return make_vector<1>(args);
  }
}

template <Int dim>
any make_matrix(ExprReader::Args& args) {
  auto v = zero_matrix<dim, dim>();
  OMEGA_H_CHECK(args.size() == std::size_t(dim * dim));
  for (Int i = 0; i < dim; ++i) {
    for (Int j = 0; j < dim; ++j) {
      auto& arg = args[std::size_t(i * dim + j)];
      v(i, j) = any_cast<Real>(arg);
    }
  }
  return v;
}

any make_matrix(LO size, Int dim, ExprReader::Args& args) {
  if (args.size() == 1 && args[0].type() == typeid(Real) &&
      any_cast<Real>(args[0]) == 0.0) {
    return Reals(size * dim * dim, 0.0);
  }
  if (args.size() != std::size_t(square(dim))) {
    throw ParserFail("Wrong number of arguments to matrix()\n");
  }
  bool has_arrays = false;
  for (auto& arg : args)
    if (arg.type() == typeid(Reals)) has_arrays = true;
  if (has_arrays) {
    std::vector<Read<Real>> arrays;
    for (Int i = 0; i < Int(args.size()); ++i) {
      auto& arg = args[std::size_t(i)];
      promote(size, dim, arg);
      arrays.push_back(any_cast<Reals>(arg));
    }
    OMEGA_H_CHECK(Int(arrays.size()) == square(dim));
    return Reals(interleave(arrays));
  } else {
    if (dim == 3) return make_matrix<3>(args);
    if (dim == 2) return make_matrix<2>(args);
    return make_matrix<1>(args);
  }
}

any make_symm(LO size, Int dim, ExprReader::Args& args) {
  if (args.size() != 1) {
    throw ParserFail("Wrong number of arguments to symm()\n");
  }
  auto arg = args[0];
  if (arg.type() != typeid(Reals)) {
    promote(size, dim, arg);
  }
  auto const in = any_cast<Reals>(arg);
  if (in.size() != size * dim * dim) {
    throw ParserFail("Argument to symm() was not sized as full tensors\n");
  }
  return matrices_to_symms(in, dim);
}

any eval_exp(LO size, ExprReader::Args& args) {
  if (args.size() != 1) {
    std::stringstream ss;
    ss << "exp() takes exactly one argument, given " << args.size();
    throw ParserFail(ss.str());
  }
  auto& in_any = args.at(0);
  if (in_any.type() == typeid(Real)) {
    return std::exp(any_cast<Real>(in_any));
  } else if (in_any.type() == typeid(Reals)) {
    auto a = any_cast<Reals>(in_any);
    if (a.size() != size) {
      throw ParserFail("exp() given array that wasn't scalars");
    }
    auto out = Write<Real>(a.size());
    auto f = OMEGA_H_LAMBDA(LO i) { out[i] = std::exp(a[i]); };
    parallel_for(a.size(), f, "eval_exp(Reals)");
    return Reals(out);
  } else {
    std::stringstream ss;
    ss << "unexpected argument type " << in_any.type().name() << " to exp()\n";
    throw ParserFail(ss.str());
  }
}

any eval_sqrt(LO size, ExprReader::Args& args) {
  if (args.size() != 1) {
    std::stringstream ss;
    ss << "sqrt() takes exactly one argument, given " << args.size();
    throw ParserFail(ss.str());
  }
  auto& in_any = args.at(0);
  if (in_any.type() == typeid(Real)) {
    return std::sqrt(any_cast<Real>(in_any));
  } else if (in_any.type() == typeid(Reals)) {
    auto a = any_cast<Reals>(in_any);
    if (a.size() != size) {
      throw ParserFail("sqrt() given array that wasn't scalars");
    }
    auto out = Write<Real>(a.size());
    auto f = OMEGA_H_LAMBDA(LO i) { out[i] = std::sqrt(a[i]); };
    parallel_for(a.size(), f, "eval_sqrt(Reals)");
    return Reals(out);
  } else {
    std::stringstream ss;
    ss << "unexpected argument type " << in_any.type().name() << " to sqrt()\n";
    throw ParserFail(ss.str());
  }
}

any eval_sin(LO size, ExprReader::Args& args) {
  if (args.size() != 1) {
    std::stringstream ss;
    ss << "sin() takes exactly one argument, given " << args.size();
    throw ParserFail(ss.str());
  }
  auto& in_any = args.at(0);
  if (in_any.type() == typeid(Real)) {
    return std::sin(any_cast<Real>(in_any));
  } else if (in_any.type() == typeid(Reals)) {
    auto a = any_cast<Reals>(in_any);
    if (a.size() != size) {
      throw ParserFail("sin() given array that wasn't scalars");
    }
    auto out = Write<Real>(a.size());
    auto f = OMEGA_H_LAMBDA(LO i) { out[i] = std::sin(a[i]); };
    parallel_for(a.size(), f, "eval_sin(Reals)");
    return Reals(out);
  } else {
    std::stringstream ss;
    ss << "unexpected argument type " << in_any.type().name() << " to sin()\n";
    throw ParserFail(ss.str());
  }
}

any eval_cos(LO size, ExprReader::Args& args) {
  if (args.size() != 1) {
    std::stringstream ss;
    ss << "cos() takes exactly one argument, given " << args.size();
    throw ParserFail(ss.str());
  }
  auto& in_any = args.at(0);
  if (in_any.type() == typeid(Real)) {
    return std::cos(any_cast<Real>(in_any));
  } else if (in_any.type() == typeid(Reals)) {
    auto a = any_cast<Reals>(in_any);
    if (a.size() != size) {
      throw ParserFail("cos() given array that wasn't scalars");
    }
    auto out = Write<Real>(a.size());
    auto f = OMEGA_H_LAMBDA(LO i) { out[i] = std::cos(a[i]); };
    parallel_for(a.size(), f, "eval_cos(Reals)");
    return Reals(out);
  } else {
    std::stringstream ss;
    ss << "unexpected argument type " << in_any.type().name() << " to cos()\n";
    throw ParserFail(ss.str());
  }
}

any eval_erf(LO size, ExprReader::Args& args) {
  if (args.size() != 1) {
    std::stringstream ss;
    ss << "erf() takes exactly one argument, given " << args.size();
    throw ParserFail(ss.str());
  }
  auto& in_any = args.at(0);
  if (in_any.type() == typeid(Real)) {
    return std::erf(any_cast<Real>(in_any));
  } else if (in_any.type() == typeid(Reals)) {
    auto a = any_cast<Reals>(in_any);
    if (a.size() != size) {
      throw ParserFail("erf() given array that wasn't scalars");
    }
    auto out = Write<Real>(a.size());
    auto f = OMEGA_H_LAMBDA(LO i) { out[i] = std::erf(a[i]); };
    parallel_for(a.size(), f, "eval_erf(Reals)");
    return Reals(out);
  } else {
    std::stringstream ss;
    ss << "unexpected argument type " << in_any.type().name() << " to erf()\n";
    throw ParserFail(ss.str());
  }
}

template <Int dim>
any eval_norm(LO size, ExprReader::Args& args) {
  if (args.size() != 1) {
    std::stringstream ss;
    ss << "norm() takes exactly one argument, given " << args.size();
    throw ParserFail(ss.str());
  }
  auto& in_any = args.at(0);
  if (in_any.type() == typeid(Vector<dim>)) {
    return Omega_h::norm(any_cast<Vector<dim>>(in_any));
  } else if (in_any.type() == typeid(Reals)) {
    auto a = any_cast<Reals>(in_any);
    if (a.size() != size * dim) {
      throw ParserFail("norm() given array that wasn't vectors");
    }
    auto out = Write<Real>(size);
    auto f = OMEGA_H_LAMBDA(LO i) {
      out[i] = Omega_h::norm(Omega_h::get_vector<dim>(a, i));
    };
    parallel_for(size, f, "eval_norm(Reals)");
    return Reals(out);
  } else {
    std::stringstream ss;
    ss << "unexpected argument type " << in_any.type().name() << " to norm()\n";
    throw ParserFail(ss.str());
  }
}

any eval_norm(Int dim, LO size, ExprReader::Args& args) {
  if (dim == 3) return eval_norm<3>(size, args);
  if (dim == 2) return eval_norm<2>(size, args);
  return eval_norm<1>(size, args);
}

}  // end anonymous namespace

ExprEnv::ExprEnv(LO size_in, Int dim_in) : size(size_in), dim(dim_in) {
  auto local_size = size;
  auto local_dim = dim;
  register_function(
      "exp", [=](Args& args) { return eval_exp(local_size, args); });
  register_function(
      "sqrt", [=](Args& args) { return eval_sqrt(local_size, args); });
  register_function(
      "sin", [=](Args& args) { return eval_sin(local_size, args); });
  register_function(
      "cos", [=](Args& args) { return eval_cos(local_size, args); });
  register_function(
      "erf", [=](Args& args) { return eval_erf(local_size, args); });
  register_function("vector",
      [=](Args& args) { return make_vector(local_size, local_dim, args); });
  register_function("symm",
      [=](Args& args) { return make_symm(local_size, local_dim, args); });
  register_function("matrix",
      [=](Args& args) { return make_matrix(local_size, local_dim, args); });
  register_function("tensor",
      [=](Args& args) { return make_matrix(local_size, local_dim, args); });
  register_function("norm",
      [=](Args& args) { return eval_norm(local_dim, local_size, args); });
  register_variable("d", any(Real(dim)));
  if (dim == 3) register_variable("I", any(identity_matrix<3, 3>()));
  if (dim == 2) register_variable("I", any(identity_matrix<2, 2>()));
  if (dim == 1) register_variable("I", any(identity_matrix<1, 1>()));
  register_variable("pi", any(Real(Omega_h::PI)));
}

void ExprEnv::register_variable(std::string const& name, any const& value) {
  OMEGA_H_TIME_FUNCTION;
  // OMEGA_H_CHECK(variables.find(name) == variables.end());
  OMEGA_H_CHECK(functions.find(name) == functions.end());
  variables[name] = value;
}

void ExprEnv::register_function(
    std::string const& name, Function const& value) {
  OMEGA_H_CHECK(variables.find(name) == variables.end());
  // OMEGA_H_CHECK(functions.find(name) == functions.end());
  functions[name] = value;
}

void ExprEnv::repeat(any& x) { promote(size, dim, x); }

ExprReader::ExprReader(LO size_in, Int dim_in)
    : Reader(math_lang::ask_reader_tables()), env(size_in, dim_in) {}

ExprReader::~ExprReader() {}

void ExprReader::register_variable(std::string const& name, any const& value) {
  env.register_variable(name, value);
}

void ExprReader::register_function(
    std::string const& name, Function const& value) {
  env.register_function(name, value);
}

void ExprReader::repeat(any& x) { env.repeat(x); }

any ExprReader::at_shift(int token, std::string& text) {
  switch (token) {
    case math_lang::TOK_NAME: {
      return text;
    }
    case math_lang::TOK_CONST: {
      return std::atof(text.c_str());
    }
  }
  return any();
}

any ExprReader::at_reduce(int prod, std::vector<any>& rhs) {
  OMEGA_H_TIME_FUNCTION;
  switch (prod) {
    case math_lang::PROD_PROGRAM: {
      if (rhs.at(1).empty()) {
        throw ParserFail(
            "Omega_h::ExprReader needs an expression to evaluate!");
      }
      return rhs.at(1);
    }
    case math_lang::PROD_NO_STATEMENTS:
    case math_lang::PROD_NO_EXPR:
    case math_lang::PROD_NEXT_STATEMENT: {
      return any();
    }
    case math_lang::PROD_ASSIGN: {
      std::string const& name = any_cast<std::string&>(rhs.at(0));
      env.variables[name] = std::move(rhs.at(4));
      return any();
    }
    case math_lang::PROD_YES_EXPR:
    case math_lang::PROD_EXPR:
    case math_lang::PROD_TERNARY_DECAY:
    case math_lang::PROD_OR_DECAY:
    case math_lang::PROD_AND_DECAY:
    case math_lang::PROD_ADD_SUB_DECAY:
    case math_lang::PROD_MUL_DIV_DECAY:
    case math_lang::PROD_POW_DECAY:
    case math_lang::PROD_NEG_DECAY:
    case math_lang::PROD_SOME_ARGS:
    case math_lang::PROD_CONST:
      return rhs.at(0);
    case math_lang::PROD_TERNARY:
      promote(env.size, env.dim, rhs.at(3), rhs.at(6));
      return ternary(env.size, env.dim, rhs.at(0), rhs.at(3), rhs.at(6));
    case math_lang::PROD_OR:
      promote(env.size, env.dim, rhs.at(0), rhs.at(3));
      return eval_or(rhs.at(0), rhs.at(3));
    case math_lang::PROD_AND:
      promote(env.size, env.dim, rhs.at(0), rhs.at(3));
      return eval_and(rhs.at(0), rhs.at(3));
    case math_lang::PROD_GT:
      promote(env.size, env.dim, rhs.at(0), rhs.at(3));
      return gt(rhs.at(0), rhs.at(3));
    case math_lang::PROD_LT:
      promote(env.size, env.dim, rhs.at(0), rhs.at(3));
      return lt(rhs.at(0), rhs.at(3));
    case math_lang::PROD_GEQ:
    case math_lang::PROD_LEQ:
      throw ParserFail("Operators <= and >= not supported yet");
    case math_lang::PROD_EQ:
      promote(env.size, env.dim, rhs.at(0), rhs.at(3));
      return eq(rhs.at(0), rhs.at(3));
    case math_lang::PROD_ADD:
      promote(env.size, env.dim, rhs.at(0), rhs.at(3));
      return add(env.dim, rhs.at(0), rhs.at(3));
    case math_lang::PROD_SUB:
      promote(env.size, env.dim, rhs.at(0), rhs.at(3));
      return sub(env.dim, rhs.at(0), rhs.at(3));
    case math_lang::PROD_MUL:
      promote(env.size, env.dim, rhs.at(0), rhs.at(3));
      return mul(env.size, env.dim, rhs.at(0), rhs.at(3));
    case math_lang::PROD_DIV:
      promote(env.size, env.dim, rhs.at(0), rhs.at(3));
      return div(env.dim, rhs.at(0), rhs.at(3));
    case math_lang::PROD_POW:
      promote(env.size, env.dim, rhs.at(0), rhs.at(3));
      return eval_pow(env.dim, rhs.at(0), rhs.at(3));
    case math_lang::PROD_CALL: {
      ScopedTimer("call");
      auto& name = any_cast<std::string&>(rhs.at(0));
      auto& args = any_cast<Args&>(rhs.at(4));
      auto vit = env.variables.find(name);
      if (vit == env.variables.end()) {
        /* function call */
        auto fit = env.functions.find(name);
        if (fit == env.functions.end()) {
          std::stringstream ss;
          ss << "\"" << name
             << "\" is neither a variable nor a function name\n";
          throw ParserFail(ss.str());
        }
        return fit->second(args);
      } else {
        /* access operator for vector/matrix */
        auto& val = vit->second;
        return access(env.size, env.dim, val, args);
      }
    }
    case math_lang::PROD_NO_ARGS: {
      return Args{};
    }
    case math_lang::PROD_FIRST_ARG: {
      Args args;
      args.push_back(std::move(rhs.at(0)));
      return any(std::move(args));
    }
    case math_lang::PROD_NEXT_ARG: {
      auto& args = any_cast<Args&>(rhs.at(0));
      args.push_back(std::move(rhs.at(3)));
      return any(std::move(args));
    }
    case math_lang::PROD_NEG:
      return neg(env.dim, rhs.at(2));
    case math_lang::PROD_VAL_PARENS:
    case math_lang::PROD_BOOL_PARENS:
      return any(std::move(rhs.at(2)));
    case math_lang::PROD_VAR: {
      auto& name = any_cast<std::string&>(rhs.at(0));
      auto it = env.variables.find(name);
      if (it == env.variables.end()) {
        std::stringstream ss;
        ss << "unknown variable name \"" << name << "\"\n";
        throw ParserFail(ss.str());
      }
      return it->second;
    }
  }
  return any();
}

ExprOp::~ExprOp() {}

struct ConstOp final : public ExprOp {
  double value;
  OpPtr rhs;
  virtual ~ConstOp() override final = default;
  ConstOp(double value_in) : value(value_in) {}
  virtual any eval(ExprEnv& env) override final;
};
any ConstOp::eval(ExprEnv&) { return value; }

struct SemicolonOp final : public ExprOp {
  OpPtr lhs;
  OpPtr rhs;
  virtual ~SemicolonOp() override final = default;
  SemicolonOp(OpPtr lhs_in, OpPtr rhs_in) : lhs(lhs_in), rhs(rhs_in) {}
  virtual any eval(ExprEnv& env) override final;
};
any SemicolonOp::eval(ExprEnv& env) {
  lhs->eval(env);  // LHS result ignored
  return rhs->eval(env);
}

struct AssignOp final : public ExprOp {
  std::string name;
  OpPtr rhs;
  virtual ~AssignOp() override final = default;
  AssignOp(std::string const& name_in, OpPtr rhs_in)
      : name(name_in), rhs(rhs_in) {}
  virtual any eval(ExprEnv& env) override final;
};
any AssignOp::eval(ExprEnv& env) {
  env.variables[name] = rhs->eval(env);
  return any();
}

struct VarOp final : public ExprOp {
  std::string name;
  virtual ~VarOp() override final = default;
  VarOp(std::string const& name_in) : name(name_in) {}
  virtual any eval(ExprEnv& env) override final;
};
any VarOp::eval(ExprEnv& env) {
  auto it = env.variables.find(name);
  if (it == env.variables.end()) {
    std::stringstream ss;
    ss << "unknown variable name \"" << name << "\"\n";
    throw ParserFail(ss.str());
  }
  return it->second;
}

struct NegOp final : public ExprOp {
  OpPtr rhs;
  virtual ~NegOp() override final = default;
  NegOp(OpPtr rhs_in) : rhs(rhs_in) {}
  virtual any eval(ExprEnv& env) override final;
};
any NegOp::eval(ExprEnv& env) { return neg(env.dim, rhs->eval(env)); }

struct TernaryOp final : public ExprOp {
  OpPtr cond;
  OpPtr lhs;
  OpPtr rhs;
  virtual ~TernaryOp() override final = default;
  TernaryOp(OpPtr cond_in, OpPtr lhs_in, OpPtr rhs_in)
      : cond(cond_in), lhs(lhs_in), rhs(rhs_in) {}
  virtual any eval(ExprEnv& env) override final;
};
any TernaryOp::eval(ExprEnv& env) {
  auto lhs_val = lhs->eval(env);
  auto rhs_val = rhs->eval(env);
  promote(env.size, env.dim, lhs_val, rhs_val);
  return ternary(env.size, env.dim, cond->eval(env), lhs_val, rhs_val);
}

struct CallOp final : public ExprOp {
  std::string name;
  std::vector<OpPtr> rhs;
  ExprEnv::Args args;
  virtual ~CallOp() override final = default;
  CallOp(std::string const& name_in, ExprEnv::Args const& args_in)
      : name(name_in) {
    for (auto& arg : args_in) {
      OpPtr op = any_cast<OpPtr>(arg);
      rhs.push_back(op);
    }
    args.reserve(rhs.size());
  }
  virtual any eval(ExprEnv& env) override final;
};
any CallOp::eval(ExprEnv& env) {
  args.resize(rhs.size());
  for (std::size_t i = 0; i < rhs.size(); ++i) {
    args[i] = rhs[i]->eval(env);
  }
  auto vit = env.variables.find(name);
  if (vit == env.variables.end()) {
    /* function call */
    auto fit = env.functions.find(name);
    if (fit == env.functions.end()) {
      std::stringstream ss;
      ss << "\"" << name << "\" is neither a variable nor a function name\n";
      throw ParserFail(ss.str());
    }
    auto result = fit->second(args);
    args.clear();
    return result;
  } else {
    /* access operator for vector/matrix */
    auto& val = vit->second;
    auto result = access(env.size, env.dim, val, args);
    args.clear();
    return result;
  }
}

#define OMEGA_H_BINARY_OP(ClassName, func_call)                                \
  struct ClassName final : public ExprOp {                                     \
    OpPtr lhs;                                                                 \
    OpPtr rhs;                                                                 \
    virtual ~ClassName() override final = default;                             \
    ClassName(OpPtr lhs_in, OpPtr rhs_in) : lhs(lhs_in), rhs(rhs_in) {}        \
    virtual any eval(ExprEnv& env) override final;                             \
  };                                                                           \
  any ClassName::eval(ExprEnv& env) {                                          \
    auto lhs_val = lhs->eval(env);                                             \
    auto rhs_val = rhs->eval(env);                                             \
    promote(env.size, env.dim, lhs_val, rhs_val);                              \
    return func_call;                                                          \
  }

OMEGA_H_BINARY_OP(OrOp, eval_or(lhs_val, rhs_val));
OMEGA_H_BINARY_OP(AndOp, eval_and(lhs_val, rhs_val));
OMEGA_H_BINARY_OP(GtOp, gt(lhs_val, rhs_val));
OMEGA_H_BINARY_OP(LtOp, lt(lhs_val, rhs_val));
OMEGA_H_BINARY_OP(EqOp, eq(lhs_val, rhs_val));
OMEGA_H_BINARY_OP(AddOp, add(env.dim, lhs_val, rhs_val));
OMEGA_H_BINARY_OP(SubOp, sub(env.dim, lhs_val, rhs_val));
OMEGA_H_BINARY_OP(MulOp, mul(env.size, env.dim, lhs_val, rhs_val));
OMEGA_H_BINARY_OP(DivOp, div(env.dim, lhs_val, rhs_val));
OMEGA_H_BINARY_OP(PowOp, eval_pow(env.dim, lhs_val, rhs_val));

#undef OMEGA_H_BINARY_OP

#define OMEGA_H_BINARY_REDUCE(ClassName)                                       \
  {                                                                            \
    OpPtr lhs_op = any_cast<OpPtr>(rhs.at(0));                                 \
    OpPtr rhs_op = any_cast<OpPtr>(rhs.at(3));                                 \
    return OpPtr(new ClassName(lhs_op, rhs_op));                               \
  }

ExprOpsReader::ExprOpsReader() : Reader(math_lang::ask_reader_tables()) {}

OpPtr ExprOpsReader::read_ops(std::string const& expr) {
  return any_cast<OpPtr>(read_string(expr, expr));
}

any ExprOpsReader::at_shift(int token, std::string& text) {
  switch (token) {
    case math_lang::TOK_NAME: {
      return text;
    }
    case math_lang::TOK_CONST: {
      return OpPtr(new ConstOp(std::atof(text.c_str())));
    }
  }
  return any();
}

any ExprOpsReader::at_reduce(int prod, std::vector<any>& rhs) {
  OMEGA_H_TIME_FUNCTION;
  using std::swap;
  switch (prod) {
    case math_lang::PROD_PROGRAM: {
      if (rhs.at(1).empty()) {
        throw ParserFail(
            "Omega_h::ExprReader needs an expression to evaluate!");
      }
      if (rhs.at(0).empty()) {
        return any(std::move(rhs.at(1)));
      } else {
        auto op_lhs = any_cast<OpPtr>(rhs.at(0));
        auto op_rhs = any_cast<OpPtr>(rhs.at(1));
        return OpPtr(new SemicolonOp(op_lhs, op_rhs));
      }
    }
    case math_lang::PROD_NO_STATEMENTS:
    case math_lang::PROD_NO_EXPR: {
      return any();
    }
    case math_lang::PROD_NEXT_STATEMENT: {
      if (rhs.at(0).empty()) {
        return any(std::move(rhs.at(1)));
      } else {
        auto op_lhs = any_cast<OpPtr>(rhs.at(0));
        auto op_rhs = any_cast<OpPtr>(rhs.at(1));
        return OpPtr(new SemicolonOp(op_lhs, op_rhs));
      }
    }
    case math_lang::PROD_ASSIGN: {
      std::string const& name = any_cast<std::string&>(rhs.at(0));
      OpPtr rhs_op = any_cast<OpPtr>(rhs.at(4));
      return OpPtr(new AssignOp(name, rhs_op));
    }
    case math_lang::PROD_YES_EXPR:
    case math_lang::PROD_EXPR:
    case math_lang::PROD_TERNARY_DECAY:
    case math_lang::PROD_OR_DECAY:
    case math_lang::PROD_AND_DECAY:
    case math_lang::PROD_ADD_SUB_DECAY:
    case math_lang::PROD_MUL_DIV_DECAY:
    case math_lang::PROD_POW_DECAY:
    case math_lang::PROD_NEG_DECAY:
    case math_lang::PROD_SOME_ARGS:
    case math_lang::PROD_CONST:
      return any(std::move(rhs.at(0)));
    case math_lang::PROD_TERNARY: {
      OpPtr cond_op = any_cast<OpPtr>(rhs.at(0));
      OpPtr lhs_op = any_cast<OpPtr>(rhs.at(3));
      OpPtr rhs_op = any_cast<OpPtr>(rhs.at(6));
      return OpPtr(new TernaryOp(cond_op, lhs_op, rhs_op));
    }
    case math_lang::PROD_OR:
      OMEGA_H_BINARY_REDUCE(OrOp)
    case math_lang::PROD_AND:
      OMEGA_H_BINARY_REDUCE(AndOp)
    case math_lang::PROD_GT:
      OMEGA_H_BINARY_REDUCE(GtOp)
    case math_lang::PROD_LT:
      OMEGA_H_BINARY_REDUCE(LtOp)
    case math_lang::PROD_GEQ:
    case math_lang::PROD_LEQ:
      throw ParserFail("Operators <= and >= not supported yet");
    case math_lang::PROD_EQ:
      OMEGA_H_BINARY_REDUCE(EqOp)
    case math_lang::PROD_ADD:
      OMEGA_H_BINARY_REDUCE(AddOp)
    case math_lang::PROD_SUB:
      OMEGA_H_BINARY_REDUCE(SubOp)
    case math_lang::PROD_MUL:
      OMEGA_H_BINARY_REDUCE(MulOp)
    case math_lang::PROD_DIV:
      OMEGA_H_BINARY_REDUCE(DivOp)
    case math_lang::PROD_POW:
      OMEGA_H_BINARY_REDUCE(PowOp)
    case math_lang::PROD_CALL: {
      auto& name = any_cast<std::string&>(rhs.at(0));
      auto& args = any_cast<ExprEnv::Args&>(rhs.at(4));
      return OpPtr(new CallOp(name, args));
    }
    case math_lang::PROD_NO_ARGS:
      return ExprEnv::Args{};
    case math_lang::PROD_FIRST_ARG: {
      ExprEnv::Args args;
      args.push_back(std::move(rhs.at(0)));
      return any(std::move(args));
    }
    case math_lang::PROD_NEXT_ARG: {
      auto& args = any_cast<ExprEnv::Args&>(rhs.at(0));
      args.push_back(std::move(rhs.at(3)));
      return any(std::move(args));
    }
    case math_lang::PROD_NEG: {
      OpPtr rhs_op = any_cast<OpPtr>(rhs.at(2));
      return OpPtr(new NegOp(rhs_op));
    }
    case math_lang::PROD_VAL_PARENS:
    case math_lang::PROD_BOOL_PARENS:
      return any(std::move(rhs.at(2)));
    case math_lang::PROD_VAR: {
      auto& name = any_cast<std::string&>(rhs.at(0));
      return OpPtr(new VarOp(name));
    }
  }
  return any();
}

std::string ExprEnv::string(int verbose) {
  //std::map<std::string, any> variables;
  //std::map<std::string, Function> functions;
  //LO size;
  //Int dim;
  std::ostringstream oss, vs;
  std::string sep = "";
  for(auto i : variables) {
    auto v = i.second;
    std::string str = any_cast<std::string>(v);
    vs << sep << "{" << i.first << " : " << str << "}";
    sep = " ";
  }
  oss << "ExprEnv:"
      << "\n  size      = " << size
      << "\n  dim       = " << dim
      << "\n  variables = " << vs.str();
  (void)verbose;
  return oss.str();
}

#undef OMEGA_H_BINARY_REDUCE

}  // end namespace Omega_h
