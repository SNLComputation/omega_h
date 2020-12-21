#ifndef OMEGA_H_DEFINES_HPP
#define OMEGA_H_DEFINES_HPP

#include <cstdint>

#include <Omega_h_macros.h>

enum Omega_h_Type {
  OMEGA_H_I8 = 0,
  OMEGA_H_I32 = 2,
  OMEGA_H_I64 = 3,
  OMEGA_H_F64 = 5,
  OMEGA_H_REAL = OMEGA_H_F64,
};

enum { OMEGA_H_DIMS = 4 };

enum Omega_h_EntDim {
  OMEGA_H_VERT = 0,
  OMEGA_H_EDGE = 1,
  OMEGA_H_FACE = 2,
  OMEGA_H_REGION = 3
};

enum Omega_h_Op { OMEGA_H_MIN, OMEGA_H_MAX, OMEGA_H_SUM };

enum Omega_h_Comparison { OMEGA_H_SAME, OMEGA_H_MORE, OMEGA_H_DIFF };

enum Omega_h_Transfer {
  OMEGA_H_INHERIT,
  OMEGA_H_LINEAR_INTERP,
  OMEGA_H_METRIC,
  OMEGA_H_DENSITY,
  OMEGA_H_CONSERVE,
  OMEGA_H_MOMENTUM_VELOCITY,
  OMEGA_H_POINTWISE,
};

enum Omega_h_Parting {
  OMEGA_H_ELEM_BASED,
  OMEGA_H_GHOSTED,
  OMEGA_H_VERT_BASED,
};

enum Omega_h_Source {
  OMEGA_H_CONSTANT,
  OMEGA_H_VARIATION,
  OMEGA_H_DERIVATIVE,
  OMEGA_H_GIVEN,
  OMEGA_H_IMPLIED,
  OMEGA_H_CURVATURE,
};

// controls the conversion of anisotropic metrics to isotropic ones
enum Omega_h_Isotropy {
  OMEGA_H_ANISOTROPIC,  // keep anisotropy
  OMEGA_H_ISO_LENGTH,   // use smallest length
  OMEGA_H_ISO_SIZE,     // use equivalent volume
};

/* determines whether a metric is allowed to be scaled
 * to satisfy an element count constraint
 */
enum Omega_h_Scales {
  OMEGA_H_ABSOLUTE,
  OMEGA_H_SCALES,
};

enum Omega_h_Family { OMEGA_H_SIMPLEX, OMEGA_H_HYPERCUBE };

namespace Omega_h {

typedef std::int8_t I8;
typedef std::int16_t I16;
typedef std::int32_t I32;
typedef std::int64_t I64;
typedef I8 Byte;
typedef I32 Int;
typedef I32 LO;
typedef I32 ClassId;
typedef I64 GO;
typedef double Real;

static constexpr Real PI = OMEGA_H_PI;
static constexpr Real EPSILON = OMEGA_H_EPSILON;
static constexpr Int MANTISSA_BITS = OMEGA_H_MANTISSA_BITS;
static constexpr bool be_verbose = true;
static constexpr bool dont_be_verbose = false;

enum { DIMS = OMEGA_H_DIMS };

enum : Int {
  VERT = OMEGA_H_VERT,
  EDGE = OMEGA_H_EDGE,
  FACE = OMEGA_H_FACE,
  REGION = OMEGA_H_REGION
};

enum { TOPO_TYPES = 8 };

}  // end namespace Omega_h

enum class Topo_type : int {
  vertex,
  edge,
  triangle,
  quadrilateral,
  tetrahedron,
  hexahedron,
  wedge,
  pyramid
};

#endif
