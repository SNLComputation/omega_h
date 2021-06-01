#ifndef OMEGA_H_DBG_HPP
#define OMEGA_H_DBG_HPP
#include "Omega_h_comm.hpp"
#include "Omega_h_fail.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>

extern Omega_h::Comm *DBG_COMM;
extern bool dbg_print_global;

#define ERROUTFL(a) do { std::ostringstream _oss_; _oss_ << __FILE__ << ":" << __LINE__ << " :dbg: " << a; OMEGA_H_CHECK(false); } while(0)
#define FILE_LINE() (std::string(__FILE__) + ":" + std::to_string(__LINE__))

#ifdef OMEGA_H_DBG

namespace Omega_h {
#ifdef OMEGA_H_USE_CUDA
OMEGA_H_INLINE std::string proc() { return ""; }
OMEGA_H_INLINE std::string rank() { return ""; }
OMEGA_H_INLINE std::string size() { return ""; }
#else
OMEGA_H_INLINE std::string proc() { std::ostringstream _oss_; _oss_ << "P" << (DBG_COMM ? DBG_COMM->rank() : 0) << ": "; return _oss_.str(); }
OMEGA_H_INLINE std::string rank() { std::ostringstream _oss_; _oss_ << (DBG_COMM ? DBG_COMM->rank() : 0); return _oss_.str(); }
OMEGA_H_INLINE std::string size() { std::ostringstream _oss_; _oss_ << (DBG_COMM ? DBG_COMM->size() : 0); return _oss_.str(); }
#endif
}

#  define TASK_0_cout if(DBG_COMM && (0 == DBG_COMM->rank())) std::cout

#  define TRACK0(a) do { std::ostringstream _oss_; _oss_ << __FILE__ << ":" << __LINE__ << " :dbg: " << #a << ": " << a; TASK_0_cout << _oss_.str() << std::endl; } while(0)
#  define TRACK1(a) do { std::ostringstream _oss_; _oss_ << __FILE__ << ":" << __LINE__ << " :dbg: " << a; TASK_0_cout << _oss_.str() << std::endl; } while(0)
#  define TRACK2(a) do { std::ostringstream _oss_; _oss_ << __FILE__ << ":" << __LINE__ << " :dbg: " << #a << ": " << a; \
                        TASK_0_cout << _oss_.str() << std::endl; if(DBG_COMM) DBG_COMM->barrier(); } while(0)
#  define TRACK() TRACK0("track")

#  define PCOUT(a) do {                                             \
    if (DBG_COMM) {                                                 \
      DBG_COMM->barrier();                                          \
      for (int _irank_=0; _irank_ < DBG_COMM->size(); _irank_++) {  \
        if (DBG_COMM->rank() == _irank_) {                          \
          std::ostringstream _oss_;                                 \
          _oss_ << "\nP" << _irank_ << ":dbg: " << a;                 \
          std::cout << _oss_.str() << std::flush;                   \
          std::cout.rdbuf()->pubsync();                             \
        }                                                           \
        DBG_COMM->barrier();                                        \
      }                                                             \
    }                                                               \
  } while(0)

#  define PCOUTFL(a) do {                                               \
    if (DBG_COMM) {                                                     \
      for (int _irank_=0; _irank_ < DBG_COMM->size(); _irank_++) {      \
        if (DBG_COMM->rank() == _irank_) {                              \
          std::ostringstream _oss_;                                     \
          _oss_ << "\nP" << _irank_ << ": " << __FILE__ << ":" << __LINE__ << " :dbg: " << a; \
          std::cout << _oss_.str() << std::flush;                       \
          std::cout.rdbuf()->pubsync();                                 \
        }                                                               \
        DBG_COMM->barrier();                                            \
      }                                                                 \
    }                                                                   \
  } while(0)

#  define PXOUTFL(a) do {                                               \
    if (DBG_COMM) {                                                     \
      std::cout << "\nP" << DBG_COMM->rank() << ": " << __FILE__ << ":" << __LINE__ << " :dbg: " << a; \
      std::cout << std::flush;                                          \
    }                                                                   \
  } while(0)

#  define PXOUT(a) do {                                                 \
    if (DBG_COMM) {                                                     \
      std::cout << "\nP" << DBG_COMM->rank() << ":dbg: " << a;            \
      std::cout << std::flush;                                          \
    }                                                                   \
  } while(0)


#else // OMEGA_H_DBG

namespace Omega_h {
OMEGA_H_INLINE std::string proc() { return "P0: "; }
OMEGA_H_INLINE std::string rank() { return "0"; }
OMEGA_H_INLINE std::string size() { return "1"; }
}

#  define TASK_0_cout std::cout

#  define TRACK0(a) do { } while(0)
#  define TRACK1(a) do { } while(0)
#  define TRACK() do { } while(0)

#  define TRACK2(a) do { } while(0)

#  define PCOUT(a) do { } while(0)
#  define PCOUTFL(a) do { } while(0)
#  define PXOUTFL(a) do { } while(0)
#  define PXOUT(a) do { } while(0)

#endif // OMEGA_H_DBG


namespace Omega_h {
template<class T=int, class U=int>
std::vector<T> range(T start, T stop, bool inclusive=false, U step=1) {
  std::vector<T> v;

  while(1) {
    if (inclusive && start > stop) break;
    if (!inclusive && start >= stop) break;
    v.push_back(start);
    start += step;
  }
  return v;
}

template<class T=int, class U=int>
std::vector<T> rangeInclusive(T start, T stop, U step=1) {
  return range<T, U>(start, stop, true, step);
}

template<class T=int, class U=int>
std::vector<T> rangeExclusive(T start, T stop, U step=1) {
  return range<T, U>(start, stop, false, step);
}

template<class T=int>
std::string to_string(const std::vector<T>& v) {
  std::ostringstream out; 
  std::string sep = "";
  for(auto i : v) {
    out << sep << i;
    sep = " ";
  }
  return out.str();
}

template<class T=int>
std::string to_string(const std::set<T>& v) {
  std::ostringstream out; 
  std::string sep = "";
  for(auto i : v) {
    out << sep << i;
    sep = " ";
  }
  return out.str();
}

template<class K, class V>
std::string to_string(const std::map<K,V>& v) {
  std::ostringstream out; 
  std::string sep = "";
  for(auto i : v) {
    out << sep << "{" << i.first << " : " << i.second << "}";
    sep = " ";
  }
  return out.str();
}

} // namespace Omega_h

#endif // OMEGA_H_DBG_HPP
