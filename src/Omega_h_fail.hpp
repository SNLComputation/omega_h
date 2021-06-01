#ifndef OMEGA_H_FAIL_HPP
#define OMEGA_H_FAIL_HPP

#include <Omega_h_config.h>
#ifdef OMEGA_H_ENABLE_DEMANGLED_STACKTRACE
#include <Omega_h_stacktrace.hpp>
#endif
#include <cassert>
#include <iostream>

#ifdef OMEGA_H_THROW
#include <exception>
#include <string>
#endif

namespace Omega_h {

#ifdef OMEGA_H_THROW
struct exception : public std::exception {
  exception(std::string const& msg_in);
  std::string msg;
  const char* what() const noexcept override;
};
#endif

void protect();	
#ifdef _MSC_VER
__declspec(noreturn)
#else
__attribute__((noreturn, format(printf, 1, 2)))
#endif
void fail(char const* format, ...);

}  // namespace Omega_h

#define Omega_h_fail Omega_h::fail

#if defined(OMEGA_H_USE_CUDA) && (defined(__clang__) || defined(_MSC_VER))
#  define OMEGA_H_CHECK(cond) assert(cond)
#elif defined(__CUDA_ARCH__)
#  define OMEGA_H_CHECK(cond) assert(cond)
#else
#  ifdef OMEGA_H_ENABLE_DEMANGLED_STACKTRACE
#    define OMEGA_H_CHECK(cond)                                                    \
      ((cond) ? ((void)0)                                                          \
              : Omega_h::fail(                                                     \
            "assertion %s failed at %s +%d\n%s\n", #cond, __FILE__, __LINE__, Omega_h::Stacktrace::demangled_stacktrace().c_str()))
#  else
#    define OMEGA_H_CHECK(cond)                                                    \
      ((cond) ? ((void)0)                                                          \
              : Omega_h::fail(                                                     \
                    "assertion %s failed at %s +%d\n", #cond, __FILE__, __LINE__))
#  endif
#endif

#ifndef NDEBUG
#  define OMEGA_H_CHECK_MSG(cond, a) do {                   \
    if (!(cond)) std::cout << "ERROR: " << a << std::endl;  \
    OMEGA_H_CHECK(cond) ;                                   \
  } while(0)
#  define OMEGA_H_CHECK_OP(l,op,r) do {                                 \
    if (!((l) op (r))) std::cout << "ERROR: " << "!(" << #l << "  " << #op << " " << #r << ") : " << "! (" << (l) << " " << #op << " " << (r) << ")" << std::endl; \
    OMEGA_H_CHECK((l) op (r));                                          \
  } while(0)
#else // NDEBUG
#  define OMEGA_H_CHECK_MSG(cond,a) OMEGA_H_CHECK(cond)
#  define OMEGA_H_CHECK_OP(l,op,r) OMEGA_H_CHECK((l) op (r))
#endif // NDEBUG

#if defined(__clang__) && !defined(OMEGA_H_USE_CUDA)
#define OMEGA_H_NORETURN(x) OMEGA_H_CHECK(false)
#else
#define OMEGA_H_NORETURN(x)                                                    \
  do {                                                                         \
    OMEGA_H_CHECK(false);                                                      \
    return x;                                                                  \
  } while (false)
#endif

#endif
