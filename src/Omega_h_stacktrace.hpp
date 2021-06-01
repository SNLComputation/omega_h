#ifndef OMEGA_H_STACKTRACE_HPP
#define OMEGA_H_STACKTRACE_HPP
#include "Omega_h_config.h"
#include <sstream>

namespace Omega_h {

#ifdef OMEGA_H_DBG
#define DO_STACKTRACE_POS 1
#else
#define DO_STACKTRACE_POS 0
#endif

// for more efficiency, less info, use the following:
#define DO_STACKTRACE_POS_USE_INT 0

#if DO_STACKTRACE_POS
#  if DO_STACKTRACE_POS_USE_INT
#    define STACKTRACE_POS_I(i) do { Stacktrace::s_position = (i) + __LINE__; } while(0)
#    define STACKTRACE_POS() STACKTRACE_POS_I(0)
#  else
#    define STACKTRACE_POS() do { std::ostringstream str; str << __FILE__ << ":" << __LINE__; Stacktrace::s_position = str.str(); } while(0)
#  endif

#else
#  define STACKTRACE_POS() do { } while(0)
#endif

#ifdef OMEGA_H_DBG
class Stacktrace {
public:
#if DO_STACKTRACE_POS
#if DO_STACKTRACE_POS_USE_INT
  static int s_position;
#else
  static std::string s_position;
#endif
#endif
  static void print_stacktrace(size_t sz=50, const std::string& msg="");
  static std::string stacktrace(size_t sz=50, const std::string& msg="");
  static std::string demangle(const std::string& x);
  static std::string demangled_stacktrace(size_t sz=50, bool also_mangled=false, const std::string& msg="");
};
#endif

}

#endif // OMEGA_H_STACKTRACE_HPP
