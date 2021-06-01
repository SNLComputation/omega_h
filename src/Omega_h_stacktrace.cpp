#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#if defined(__GNUC__)
#include <cxxabi.h>
#endif

#include <string>
#include <iostream>

#include "Omega_h_stacktrace.hpp"

#ifdef OMEGA_H_DBG
namespace Omega_h {

#if DO_STACKTRACE_POS
#if DO_STACKTRACE_POS_USE_INT
  int Stacktrace::s_position = 0;
#else
  std::string Stacktrace::s_position = "";
#endif
#endif

  void Stacktrace::print_stacktrace(size_t sz, const std::string& msg)
  {
    void *array[sz];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, sz);

    // print out all the frames to stderr
    fprintf(stderr, "stacktrace: message= %s:\n", msg.c_str());
    //if (!get_rank())
    backtrace_symbols_fd(array, size, 2);
  }

  std::string Stacktrace::stacktrace(size_t sz, const std::string& msg)
  {
    void *array[sz];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, sz);

    std::string ret = "stacktrace: "+msg + "\n";
    char ** syms = backtrace_symbols(array, size);
    //std::cout << "stacktrace size= " << size << std::endl;
    for (size_t i =0; i < size; ++i)
      {
        ret += std::string(syms[i])+"\n";
      }
    free(syms);
    return ret;
  }

  std::string Stacktrace::demangle(const std::string& x)
  {

    std::string ret = x;
#if defined(__GNUC__)
    std::string y = x;
    size_t lpos = y.find("(");
    size_t rpos = y.find("+");
    if (lpos != std::string::npos && rpos != std::string::npos)
      {
        y = y.substr(lpos+1,rpos-lpos-1);
      }
    //std::cout << "x= " << x << "\ny= " << y << " lpos= " << lpos << " rpos= " << rpos << " npos= " << std::string::npos << std::endl;
    int status=0;
    char *realname=0;
    realname = abi::__cxa_demangle(y.c_str(), 0, 0, &status);
    if (status != 0)
      ret = y;
    else
      {
        ret = realname;
        free(realname);
      }
#endif
    return ret;
  }

  std::string Stacktrace::demangled_stacktrace(size_t sz, bool also_mangled, const std::string& msg)
  {
    (void)msg;
    void *array[sz];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, sz);

    std::ostringstream str;
#if DO_STACKTRACE_POS
    //str << "stacktrace: pos= " << s_position << " " << msg << "\n";
#endif
    std::string ret = str.str();
    char ** syms = backtrace_symbols(array, size);
    //std::cout << "stacktrace size= " << size << std::endl;
    if (also_mangled)
      {
        ret += "\n";
        for (size_t i =0; i < size; ++i)
          {
            std::string ss = syms[i];
            ret += ss+"\n";
          }
      }
    ret += "\n";
    for (size_t i =0; i < size; ++i)
      {
        std::string ss = syms[i];
        ret += demangle(ss)+"\n";
      }
    ret += "\n";
    free(syms);
    return ret;
  }

} // namespace Omega_h
#endif // OMEGA_H_DBG
