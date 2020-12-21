#ifndef OMEGA_H_ADJ_HPP
#define OMEGA_H_ADJ_HPP

#include <Omega_h_few.hpp>
#include <Omega_h_graph.hpp>
#include <Omega_h_matrix.hpp>
#include <Omega_h_vector.hpp>

namespace Omega_h {

struct Adj : public Graph {
  OMEGA_H_INLINE Adj() {}
  explicit Adj(LOs ab2b_) : Graph(ab2b_) {}
  Adj(LOs ab2b_, Read<I8> codes_) : Graph(ab2b_), codes(codes_) {}
  Adj(LOs a2ab_, LOs ab2b_, Read<I8> codes_)
      : Graph(a2ab_, ab2b_), codes(codes_) {}
  Adj(LOs a2ab_, LOs ab2b_) : Graph(a2ab_, ab2b_) {}
  Adj(Graph g) : Graph(g) {}
  Read<I8> codes;
};

struct Parents {
  OMEGA_H_INLINE Parents() {}
  Parents(LOs parent_idx_, Read<I8> codes_)
      : parent_idx(parent_idx_), codes(codes_) {}
  LOs parent_idx;
  Read<I8> codes;
};

struct Children : public Graph {
  OMEGA_H_INLINE Children() {}
  Children(LOs a2ab_, LOs ab2b_, Read<I8> codes_)
      : Graph(a2ab_, ab2b_), codes(codes_) {}
  Read<I8> codes;
};

void find_matches(Omega_h_Family const family, Int const dim, LOs const av2v,
    LOs const bv2v, Adj const v2b, Write<LO>* a2b_out, Write<I8>* codes_out);

void find_matches(Topo_type const ent_type, LOs const av2v,
    LOs const bv2v, Adj const v2b, Write<LO>* a2b_out, Write<I8>* codes_out);

Adj reflect_down(LOs const hv2v, LOs const lv2v, Adj const v2l,
    Omega_h_Family const family, Int const high_dim, Int const low_dim);

Adj reflect_down(LOs const hv2v, LOs const lv2v, Adj const v2l,
    Topo_type const high_type, Topo_type const low_type);

Adj unmap_adjacency(LOs const a2b, Adj const b2c);

/* Given a downward adjacency, derive its corresponding upward adjacency.
   The list of upward adjacent entities will be sorted by the local
   index of the upward adjacent entity */
Adj invert_adj(Adj const down, Int const nlows_per_high, LO const nlows,
    Int const high_dim, Int const low_dim);

Adj invert_adj(Adj const down, Int const nlows_per_high, LO const nlows,
    Topo_type high_type, Topo_type low_type);

Children invert_parents(Parents const children2parents, Int const parent_dim,
    Int const nparent_dim_ents);

/* given the vertex lists for high entities,
   create vertex lists for all uses of low
   entities by high entities */
LOs form_uses(LOs const hv2v, Omega_h_Family const family, Int const high_dim,
    Int const low_dim);

LOs form_uses(LOs const hv2v, Topo_type const high_type, Topo_type const low_type);

LOs find_unique(LOs const hv2v, Omega_h_Family const family, Int const high_dim,
    Int const low_dim);

/* find_unique if not used for mixed till now because
    multiple down adjs from and below face2edge will be created.
    If found a way to efficiently combine these, build from elems2verts
    can be written
*/
LOs find_unique(LOs const hv2v, Topo_type const high_type, Topo_type const low_type);

/* for each entity (or entity use), sort its vertex list
   and express the sorting transformation as an alignment code */
template <typename T>
Read<I8> get_codes_to_canonical(Int const deg, Read<T> const ev2v);

Read<I8> find_canonical_jumps(
    Int const deg, LOs const canon, LOs const e_sorted2e);

/* given entity uses and unique entities,
   both defined by vertex lists, match
   uses to unique entities and derive their
   respective alignment codes.

   even though this is a downward adjacency, we'll
   define the code as describing how to transform
   the boundary entity into the entity use,
   since typically data is being pulled into an element
   from its boundary
*/
template <typename T>
void find_matches_ex(Int const deg, LOs const a2fv, Read<T> const av2v,
    Read<T> const bv2v, Adj const v2b, Write<LO>* a2b_out, Write<I8>* codes_out,
    bool const allow_duplicates = false);

/* for testing only, internally computes upward
   adjacency */
Adj reflect_down(LOs const hv2v, LOs const lv2v, Omega_h_Family const family,
    LO const nv, Int const high_dim, Int const low_dim);

Adj transit(Adj const h2m, Adj const m2l, Omega_h_Family const family,
    Int const high_dim, Int const low_dim);

Adj transit(Adj const h2m, Adj const m2l, Topo_type const high_type,
    Topo_type const low_type, Topo_type const mid_type);

Graph verts_across_edges(Adj const e2v, Adj const v2e);
Graph edges_across_tris(Adj const f2e, Adj const e2f);
Graph edges_across_tets(Adj const r2e, Adj const e2r);
Graph elements_across_sides(Int const dim, Adj const elems2sides,
    Adj const sides2elems, Read<I8> const side_is_exposed);

template <Int nhhl>
OMEGA_H_DEVICE Few<LO, nhhl> gather_down(LOs const& hl2l, Int h) {
  Few<LO, nhhl> hhl2l;
  for (Int i = 0; i < nhhl; ++i) {
    auto hl = h * nhhl + i;
    hhl2l[i] = hl2l[hl];
  }
  return hhl2l;
}

template <Int neev>
OMEGA_H_DEVICE Few<LO, neev> gather_verts(LOs const& ev2v, Int e) {
  return gather_down<neev>(ev2v, e);
}

template <Int neev, typename T>
OMEGA_H_DEVICE Few<T, neev> gather_values(Read<T> const& a, Few<LO, neev> v) {
  Few<T, neev> x;
  for (Int i = 0; i < neev; ++i) x[i] = a[v[i]];
  return x;
}

template <Int neev>
OMEGA_H_DEVICE Vector<neev> gather_scalars(
    Read<Real> const& a, Few<LO, neev> v) {
  Vector<neev> x;
  for (Int i = 0; i < neev; ++i) x[i] = a[v[i]];
  return x;
}

template <Int neev, Int dim>
OMEGA_H_DEVICE Matrix<dim, neev> gather_vectors(
    Reals const& a, Few<LO, neev> v) {
  Matrix<dim, neev> x;
  for (Int i = 0; i < neev; ++i) x[i] = get_vector<dim>(a, v[i]);
  return x;
}

template <Int neev, Int dim>
OMEGA_H_DEVICE Few<Matrix<dim, dim>, neev> gather_symms(
    Reals const& a, Few<LO, neev> v) {
  Few<Matrix<dim, dim>, neev> x;
  for (Int i = 0; i < neev; ++i) x[i] = get_symm<dim>(a, v[i]);
  return x;
}

#define INST_DECL(T)                                                           \
  extern template Read<I8> get_codes_to_canonical(Int deg, Read<T> ev2v);      \
  extern template void find_matches_ex(Int deg, LOs a2fv, Read<T> av2v,        \
      Read<T> bv2v, Adj v2b, Write<LO>* a2b_out, Write<I8>* codes_out, bool);
INST_DECL(LO)
INST_DECL(GO)
#undef INST_DECL

}  // end namespace Omega_h

#endif
