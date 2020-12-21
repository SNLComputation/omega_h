#ifndef OMEGA_H_ELEMENT_HPP
#define OMEGA_H_ELEMENT_HPP

#include <Omega_h_config.h>
#include <Omega_h_hypercube.hpp>
#include <Omega_h_simplex.hpp>

namespace Omega_h {

OMEGA_H_INLINE Int element_down_template(Omega_h_Family family, Int elem_dim,
    Int bdry_dim, Int which_bdry, Int which_vert) {
  return (family == OMEGA_H_SIMPLEX ? simplex_down_template(elem_dim, bdry_dim,
                                          which_bdry, which_vert)
                                    : hypercube_down_template(elem_dim,
                                          bdry_dim, which_bdry, which_vert));
}

constexpr OMEGA_H_INLINE Int element_down_template(
    Int elem_type, Int bdry_type, Int which_bdry, Int which_vert) {
  // clang-format off
  return (elem_type == 7 ?
           (bdry_type == 3 ?
             (which_bdry == 0 ?
               (which_vert == 0 ? 1 :
               (which_vert == 1 ? 0 :
               (which_vert == 2 ? 3 :
               (which_vert == 3 ? 2 : -1)))) : -1) :
           (bdry_type == 2 ?
             (which_bdry == 0 ?
               (which_vert == 0 ? 0 :
               (which_vert == 1 ? 1 :
               (which_vert == 2 ? 4 : -1))) :
             (which_bdry == 1 ?
               (which_vert == 0 ? 1 :
               (which_vert == 1 ? 2 :
               (which_vert == 2 ? 4 : -1))) :
             (which_bdry == 2 ?
               (which_vert == 0 ? 2 :
               (which_vert == 1 ? 3 :
               (which_vert == 2 ? 4 : -1))) :
             (which_bdry == 3 ?
               (which_vert == 0 ? 3 :
               (which_vert == 1 ? 0 :
               (which_vert == 2 ? 4 : -1))) : -1)))) :
           (bdry_type == 1 ?
             (which_bdry == 0 ?
               (which_vert == 0 ? 0 :
               (which_vert == 1 ? 1 : -1)) :
             (which_bdry == 1 ?
               (which_vert == 0 ? 1 :
               (which_vert == 1 ? 2 : -1)) :
             (which_bdry == 2 ?
               (which_vert == 0 ? 2 :
               (which_vert == 1 ? 3 : -1)) :
             (which_bdry == 3 ?
               (which_vert == 0 ? 3 :
               (which_vert == 1 ? 0 : -1)) :
             (which_bdry == 4 ?
               (which_vert == 0 ? 0 :
               (which_vert == 1 ? 4 : -1)) :
             (which_bdry == 5 ?
               (which_vert == 0 ? 1 :
               (which_vert == 1 ? 4 : -1)) :
             (which_bdry == 6 ?
               (which_vert == 0 ? 2 :
               (which_vert == 1 ? 4 : -1)) :
             (which_bdry == 7 ?
               (which_vert == 0 ? 3 :
               (which_vert == 1 ? 4 : -1)) : -1)))))))) :
           (bdry_type == 0 ?
             (which_bdry == 0 ? 0 :
             (which_bdry == 1 ? 1 :
             (which_bdry == 2 ? 2 :
             (which_bdry == 3 ? 3 :
             (which_bdry == 4 ? 4 : -1))))) : -1)))) :
         (elem_type == 6 ?
           (bdry_type == 3 ?
             (which_bdry == 0 ?
               (which_vert == 0 ? 0 :
               (which_vert == 1 ? 1 :
               (which_vert == 2 ? 4 :
               (which_vert == 3 ? 3 : -1)))) :
             (which_bdry == 1 ?
               (which_vert == 0 ? 1 :
               (which_vert == 1 ? 2 :
               (which_vert == 2 ? 5 :
               (which_vert == 3 ? 4 : -1)))) :
             (which_bdry == 2 ?
               (which_vert == 0 ? 2 :
               (which_vert == 1 ? 0 :
               (which_vert == 2 ? 3 :
               (which_vert == 3 ? 5 : -1)))) : -1))) :
           (bdry_type == 2 ?
             (which_bdry == 0 ?
               (which_vert == 0 ? 1 :
               (which_vert == 1 ? 0 :
               (which_vert == 2 ? 2 : -1))) :
             (which_bdry == 1 ?
               (which_vert == 0 ? 3 :
               (which_vert == 1 ? 4 :
               (which_vert == 2 ? 5 : -1)))  : -1)) :
           (bdry_type == 1 ?
             (which_bdry == 0 ?
               (which_vert == 0 ? 0 :
               (which_vert == 1 ? 1 : -1)) :
             (which_bdry == 1 ?
               (which_vert == 0 ? 1 :
               (which_vert == 1 ? 2 : -1)) :
             (which_bdry == 2 ?
               (which_vert == 0 ? 2 :
               (which_vert == 1 ? 0 : -1)) :
             (which_bdry == 3 ?
               (which_vert == 0 ? 0 :
               (which_vert == 1 ? 3 : -1)) :
             (which_bdry == 4 ?
               (which_vert == 0 ? 1 :
               (which_vert == 1 ? 4 : -1)) :
             (which_bdry == 5 ?
               (which_vert == 0 ? 2 :
               (which_vert == 1 ? 5 : -1)) :
             (which_bdry == 6 ?
               (which_vert == 0 ? 4 :
               (which_vert == 1 ? 3 : -1)) :
             (which_bdry == 7 ?
               (which_vert == 0 ? 5 :
               (which_vert == 1 ? 4 : -1)) :
             (which_bdry == 8 ?
               (which_vert == 0 ? 3 :
               (which_vert == 1 ? 5 : -1)) : -1))))))))) :
           (bdry_type == 0 ?
             (which_bdry == 0 ? 0 :
             (which_bdry == 1 ? 1 :
             (which_bdry == 2 ? 2 :
             (which_bdry == 3 ? 3 :
             (which_bdry == 4 ? 4 :
             (which_bdry == 5 ? 5 : -1)))))) : -1)))) :
         (elem_type == 5 ?
           (bdry_type == 3 ?
             hypercube_down_template(3, 2, which_bdry, which_vert) :
           (bdry_type == 1 ?
             hypercube_down_template(3, 1, which_bdry, which_vert) :
           (bdry_type == 0 ?
             hypercube_down_template(3, 0, which_bdry, which_vert) : -1))) :
         (elem_type == 4 ?
           (bdry_type == 2 ?
             simplex_down_template(3, 2, which_bdry, which_vert) :
           (bdry_type == 1 ?
             simplex_down_template(3, 1, which_bdry, which_vert) :
           (bdry_type == 0 ?
             simplex_down_template(3, 0, which_bdry, which_vert) : -1))) :
         (elem_type == 3 ?
           (bdry_type == 1 ?
             hypercube_down_template(2, 1, which_bdry, which_vert) :
           (bdry_type == 0 ?
             hypercube_down_template(2, 0, which_bdry, which_vert) : -1)) :
         (elem_type == 2 ?
           (bdry_type == 1 ?
             simplex_down_template(2, 1, which_bdry, which_vert) :
           (bdry_type == 0 ?
             simplex_down_template(2, 0, which_bdry, which_vert) : -1)) :
         (elem_type == 1 ?
           simplex_down_template(1, 0, which_bdry, which_vert) : -1)))))));
  // clang-format on
}

OMEGA_H_INLINE Int element_degree(
    Omega_h_Family family, Int from_dim, Int to_dim) {
  return (family == OMEGA_H_SIMPLEX ? simplex_degree(from_dim, to_dim)
                                    : hypercube_degree(from_dim, to_dim));
}

OMEGA_H_INLINE Int element_degree(Topo_type from_type, Topo_type to_type) {
  const int from = int(from_type);
  const int to = int(to_type);
  return (from == 0 ?
	   (to == 0 ? 1 : -1) :
	 (from == 1 ?
	   (to == 0 ? 2 :
	   (to == 1 ? 1 : -1)) :
         (from == 2 ?
	   (to == 0 ? 3 :
	   (to == 1 ? 3 :
	   (to == 2 ? 1 : -1))) :
	 (from == 3 ?
	   (to == 0 ? 4 :
	   (to == 1 ? 4 :
	   (to == 3 ? 1 : -1))) :
	 (from == 4 ?
	   (to == 0 ? 4 :
	   (to == 1 ? 6 :
	   (to == 2 ? 4 :
	   (to == 4 ? 1 : -1)))) : 
	 (from == 5 ?
	   (to == 0 ? 8 :
	   (to == 1 ? 12 :
	   (to == 3 ? 6 :
	   (to == 5 ? 1 : -1)))) :
	 (from == 6 ?
	   (to == 0 ? 6 :
	   (to == 1 ? 9 :
	   (to == 2 ? 2 :
	   (to == 3 ? 3 :
	   (to == 6 ? 1 : -1))))) :
	 (from == 7 ?
	   (to == 0 ? 5 :
	   (to == 1 ? 8 :
	   (to == 2 ? 4 :
	   (to == 3 ? 1 :
	   (to == 7 ? 1 : -1))))) : -1))))))));
}

OMEGA_H_INLINE TemplateUp element_up_template(Omega_h_Family family,
    Int elem_dim, Int bdry_dim, Int which_bdry, Int which_up) {
  return (
      family == OMEGA_H_SIMPLEX
          ? simplex_up_template(elem_dim, bdry_dim, which_bdry, which_up)
          : hypercube_up_template(elem_dim, bdry_dim, which_bdry, which_up));
}

OMEGA_H_INLINE TemplateUp element_up_template(Int elem_type, Int bdry_type, Int which_bdry, Int which_up) {
  switch (elem_type) {
    case 2:
      return simplex_up_template(2, bdry_type, which_bdry, which_up);//return simplex_up_template for tri-to-vtx
    case 3:
      return hypercube_up_template(2, bdry_type, which_bdry, which_up);//return hypercube_up_temp for quad-to-vtx
    case 4:
      return simplex_up_template(3, bdry_type, which_bdry, which_up);//return simplex_up_template for tet-to-vtx/edge
    case 5:
      return hypercube_up_template(3, bdry_type, which_bdry, which_up);//return hypercube_up_temp for hex-to-vtx/edge
    case 6:
    //new template for wedge-to-vtx & wedge-to-edge via quads
    //orderings as per simmetrix
      switch (bdry_type) {
        case 0:
          switch (which_bdry) {
            case 0:
              switch (which_up) {
                case 0:
                  return {0, 0, 0}; 
              }
              return {-1, -1, 0}; 
            case 1:
              switch (which_up) {
                case 0:
                  return {1, 0, 0}; 
              }
              return {-1, -1, 0}; 
            case 2:
              switch (which_up) {
                case 0:
                  return {2, 0, 0}; 
              }
              return {-1, -1, 0}; 
            case 3:
              switch (which_up) {
                case 0:
                  return {6, 1, 0};
              }
              return {-1, -1, 0};
	    case 4:
              switch (which_up) {
                case 0:
                  return {7, 1, 0};
              }
              return {-1, -1, 0};
	    case 5:
              switch (which_up) {
                case 0:
                  return {8, 1, 0};
              }
              return {-1, -1, 0};
          }
          return {-1, -1, 0};
        case 1:
	//up-entity is quad
          switch (which_bdry) {
            case 0:
              switch (which_up) {
                case 0:
                  return {0, 0, 0}; 
              }   
              return {-1, -1, 0}; 
            case 1:
              switch (which_up) {
                case 0:
                  return {1, 0, 0}; 
              }
              return {-1, -1, 0}; 
            case 2:
              switch (which_up) {
                case 0:
                  return {2, 0, 0}; 
              }
              return {-1, -1, 0};
            case 3:
              switch (which_up) {
                case 0:
                  return {0, 3, 0};
              }
              return {-1, -1, 0};
            case 4:
              switch (which_up) {
                case 0:
                  return {1, 3, 0};
              }
              return {-1, -1, 0};
            case 5:
              switch (which_up) {
                case 0:
                  return {2, 3, 0};
              }
              return {-1, -1, 0};
            case 6:
              switch (which_up) {
                case 0:
                  return {0, 2, 0};
              }
              return {-1, -1, 0};
            case 7:
              switch (which_up) {
                case 0:
                  return {1, 2, 0};
              }
              return {-1, -1, 0};
            case 8:
              switch (which_up) {
                case 0:
                  return {2, 2, 0};
              }
              return {-1, -1, 0};
          }
          return {-1, -1, 0};
      }
      return {-1, -1, 0};
    
    case 7: // make new template for pyramid-to-vtx & pyramid-to-edge via tris
      switch (bdry_type) {
        case 0:
          switch (which_bdry) {
            case 0:
              switch (which_up) {
                case 0:
                  return {0, 0, 0}; 
              }
              return {-1, -1, 0}; 
            case 1:
              switch (which_up) {
                case 0:
                  return {1, 0, 0}; 
              }
              return {-1, -1, 0}; 
            case 2:
              switch (which_up) {
                case 0:
                  return {2, 0, 0}; 
              }
              return {-1, -1, 0}; 
            case 3:
              switch (which_up) {
                case 0:
                  return {3, 0, 0};
              }
              return {-1, -1, 0};
	    case 4:
              switch (which_up) {
                case 0:
                  return {4, 0, 0};
              }
              return {-1, -1, 0};
          }
          return {-1, -1, 0};
        case 1:
	//up-entity is tri
          switch (which_bdry) {
            case 0:
              switch (which_up) {
                case 0:
                  return {0, 0, 0}; 
              }   
              return {-1, -1, 0}; 
            case 1:
              switch (which_up) {
                case 0:
                  return {1, 0, 0}; 
              }
              return {-1, -1, 0}; 
            case 2:
              switch (which_up) {
                case 0:
                  return {2, 0, 0}; 
              }
              return {-1, -1, 0};
            case 3:
              switch (which_up) {
                case 0:
                  return {3, 0, 0};
              }
              return {-1, -1, 0};
            case 4:
              switch (which_up) {
                case 0:
                  return {0, 2, 0};
              }
              return {-1, -1, 0};
            case 5:
              switch (which_up) {
                case 0:
                  return {1, 2, 0};
              }
              return {-1, -1, 0};
            case 6:
              switch (which_up) {
                case 0:
                  return {2, 2, 0};
              }
              return {-1, -1, 0};
            case 7:
              switch (which_up) {
                case 0:
                  return {3, 2, 0};
              }
              return {-1, -1, 0};
          }
          return {-1, -1, 0};
      }
      return {-1, -1, 0};
  }
  return {-1, -1, 0};
};

constexpr char const* dimensional_singular_name(Int dim) {
  return (dim == 3 ? "region"
                   : (dim == 2 ? "face"
                               : (dim == 1 ? "edge"
                                           : (dim == 0 ? "vertex" : nullptr))));
}

constexpr char const* dimensional_singular_name(Topo_type ent_type) {
  return (int(ent_type) == 7 ? "pyramid"
	 : (int(ent_type) == 6 ? "wedge"
	   : (int(ent_type) == 5 ? "hexahedron"
	     : (int(ent_type) == 4 ? "tetrahedron"
	       : (int(ent_type) == 3 ? "quadrilateral"
		 : (int(ent_type) == 2 ? "triangle"
		   : (int(ent_type) == 1 ? "edge"
		     : (int(ent_type) == 0 ? "vertex" : nullptr))))))));
}

constexpr char const* dimensional_plural_name(Int dim) {
  return (dim == 3
              ? "regions"
              : (dim == 2 ? "faces"
                          : (dim == 1 ? "edges"
                                      : (dim == 0 ? "vertices" : nullptr))));
}

constexpr char const* dimensional_plural_name(Topo_type ent_type) {
  return (int(ent_type) == 7 ? "pyramids"
	 : (int(ent_type) == 6 ? "wedges"
	   : (int(ent_type) == 5 ? "hexahedrons"
	     : (int(ent_type) == 4 ? "tetrahedrons"
	       : (int(ent_type) == 3 ? "quadrilaterals"
		 : (int(ent_type) == 2 ? "triangles"
		   : (int(ent_type) == 1 ? "edges"
		     : (int(ent_type) == 0 ? "vertices" : nullptr))))))));
}

constexpr char const* topological_singular_name(
    Omega_h_Family family, Int dim) {
  return (family == OMEGA_H_SIMPLEX ? simplex_singular_name(dim)
                                    : hypercube_singular_name(dim));
}

constexpr char const* topological_plural_name(Omega_h_Family family, Int dim) {
  return (family == OMEGA_H_SIMPLEX ? simplex_plural_name(dim)
                                    : hypercube_plural_name(dim));
}

}  // namespace Omega_h

#endif
