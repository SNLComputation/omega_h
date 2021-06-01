#ifndef OMEGA_H_EXPR_HPP
#define OMEGA_H_EXPR_HPP

#include <functional>
#include <map>
#include <vector>

#include <Omega_h_any.hpp>
#include <Omega_h_array.hpp>
#include <Omega_h_matrix.hpp>
#include <Omega_h_reader.hpp>

namespace Omega_h {

struct ExprEnv {
  ExprEnv() = default;
  ExprEnv(LO size_in, Int dim_in);
  using Args = std::vector<any>;
  using Function = std::function<any(Args&)>;
  void register_variable(std::string const& name, any const& value);
  void register_function(std::string const& name, Function const& value);
  void repeat(any& x);
  std::map<std::string, any> variables;
  std::map<std::string, Function> functions;
  LO size;
  Int dim;
  std::string string(int verbose=0);
};

struct ExprOp {
  virtual ~ExprOp();
  virtual any eval(ExprEnv& env) = 0;
};

using OpPtr = std::shared_ptr<ExprOp>;

class ExprOpsReader final : public Reader {
 public:
  ExprOpsReader();
  virtual ~ExprOpsReader() override final = default;
  OpPtr read_ops(std::string const& str);

 protected:
  any at_shift(int token, std::string& text) override final;
  any at_reduce(int token, std::vector<any>& rhs) override final;
};

class ExprReader final : public Reader {
 public:
  using Args = ExprEnv::Args;
  using Function = ExprEnv::Function;
  ExprReader(LO size_in, Int dim_in);
  virtual ~ExprReader() override final;
  void register_variable(std::string const& name, any const& value);
  void register_function(std::string const& name, Function const& value);
  void repeat(any& x);

 protected:
  any at_shift(int token, std::string& text) override final;
  any at_reduce(int token, std::vector<any>& rhs) override final;
  ExprEnv env;
};

}  // end namespace Omega_h

#endif
