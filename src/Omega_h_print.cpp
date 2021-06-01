#include <Omega_h_print.hpp>

namespace Omega_h {

template <class T>
std::ostream& operator<<(std::ostream& stream, HostRead<T> hr) {
  stream << '{';
  using PT = promoted_t<T>;
  const bool do_designators = (hr.size() > 8);
  for (LO i = 0; i < hr.size(); ++i) {
    if (i > 0) stream << ", ";
    if (do_designators) stream << '[' << i << "]=";
    stream << PT(hr[i]);
  }
  stream << '}';
  return stream;
}

template <class T>
std::ostream& operator<<(std::ostream& stream, Read<T> r) {
  return stream << HostRead<T>(r);
}

template <class T>
std::ostream& operator<<(std::ostream& stream, Write<T> r) {
  return stream << Omega_h::read(r);
}

#define OMEGA_H_EXPL_INST(T)                                                   \
  template std::ostream& operator<<(std::ostream&, HostRead<T>);               \
  template std::ostream& operator<<(std::ostream&, Read<T>);                   \
  template std::ostream& operator<<(std::ostream&, Write<T>);
OMEGA_H_EXPL_INST(I8)
OMEGA_H_EXPL_INST(I32)
OMEGA_H_EXPL_INST(I64)
OMEGA_H_EXPL_INST(Real)
#undef OMEGA_H_EXPL_INST

}  // namespace Omega_h
