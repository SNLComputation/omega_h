#include "Omega_h_adj.hpp"

#include "Omega_h_align.hpp"
#include "Omega_h_amr.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_int_scan.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_sort.hpp"
#include "Omega_h_timer.hpp"

namespace Omega_h {

Adj unmap_adjacency(LOs const a2b, Adj const b2c) {
  OMEGA_H_TIME_FUNCTION;
  auto const b2bc = b2c.a2ab;
  auto const bc2c = b2c.ab2b;
  auto const bc_codes = b2c.codes;
  auto const b_degrees = get_degrees(b2bc);
  LOs a_degrees = unmap(a2b, b_degrees, 1);
  auto const a2ac = offset_scan(a_degrees);
  auto const na = a2b.size();
  auto const nac = a2ac.last();
  Write<LO> ac2c(nac);
  auto const ac_codes = Write<I8>(nac);
  auto f = OMEGA_H_LAMBDA(LO a) {
    auto const b = a2b[a];
    auto bc = b2bc[b];
    for (auto ac = a2ac[a]; ac < a2ac[a + 1]; ++ac) {
      ac2c[ac] = bc2c[bc];
      ac_codes[ac] = bc_codes[bc];
      ++bc;
    }
  };
  parallel_for(na, std::move(f));
  return Adj(a2ac, ac2c, ac_codes);
}

template <Int deg>
struct IsFlipped;

template <>
struct IsFlipped<4> {
  template <typename T>
  OMEGA_H_INLINE static bool is(T const adj[]) {
    return adj[3] < adj[1];
  }
};

template <>
struct IsFlipped<3> {
  template <typename T>
  OMEGA_H_INLINE static bool is(T const adj[]) {
    return adj[2] < adj[1];
  }
};

template <>
struct IsFlipped<2> {
  template <typename T>
  OMEGA_H_INLINE static bool is(T const adj[]) {
    (void)adj;
    return false;
  }
};

template <Int deg, typename T>
Read<I8> get_codes_to_canonical_deg(Read<T> const ev2v) {
  auto const nev = ev2v.size();
  auto const ne = divide_no_remainder(nev, deg);
  Write<I8> codes(ne);
  auto f = OMEGA_H_LAMBDA(LO e) {
    auto const begin = e * deg;
    /* find the smallest vertex */
    Int min_j = 0;
    auto min_v = ev2v[begin];
    for (Int j = 1; j < deg; ++j) {
      auto const ev = j + begin;
      auto const v = ev2v[ev];
      if (v < min_v) {
        min_j = j;
        min_v = v;
      }
    }
    /* rotate to make it first */
    auto const rotation = rotation_to_first(deg, min_j);
    T tmp[deg];
    rotate_adj<deg>(rotation, ev2v, begin, tmp, 0);
    auto const is_flipped = IsFlipped<deg>::is(tmp);
    codes[e] = make_code(is_flipped, rotation, 0);
  };
  parallel_for(ne, std::move(f));
  return codes;
}

template <typename T>
Read<I8> get_codes_to_canonical(Int const deg, Read<T> const ev2v) {
  OMEGA_H_TIME_FUNCTION;
  if (deg == 4) return get_codes_to_canonical_deg<4>(ev2v);
  if (deg == 3) return get_codes_to_canonical_deg<3>(ev2v);
  if (deg == 2) return get_codes_to_canonical_deg<2>(ev2v);
  OMEGA_H_NORETURN(Read<I8>());
}

/* check whether adjacent lists of (deg) vertices
   are the same */
OMEGA_H_DEVICE static bool are_equal(
    Int const deg, LOs const& canon, LO e0, LO e1) {
  auto const a = e0 * deg;
  auto const b = e1 * deg;
  for (LO j = 0; j < deg; ++j) {
    if (canon[a + j] != canon[b + j]) return false;
  }
  return true;
}

Read<I8> find_canonical_jumps(
    Int const deg, LOs const canon, LOs const e_sorted2e) {
  OMEGA_H_TIME_FUNCTION;
  auto const ne = e_sorted2e.size();
  Write<I8> jumps(ne, 0);
  auto f = OMEGA_H_LAMBDA(LO e_sorted) {
    auto const e0 = e_sorted2e[e_sorted];
    auto const e1 = e_sorted2e[e_sorted + 1];
    if (!are_equal(deg, canon, e0, e1)) jumps[e_sorted] = 1;
  };
  parallel_for(ne - 1, std::move(f));
  if (jumps.size()) jumps.set(jumps.size() - 1, 1);
  return jumps;
}

static LOs find_unique_deg(Int const deg, LOs const uv2v) {
  OMEGA_H_TIME_FUNCTION;
  auto const codes = get_codes_to_canonical(deg, uv2v);
  auto const uv2v_canon = align_ev2v(deg, uv2v, codes);
  auto const sorted2u = sort_by_keys(uv2v_canon, deg);
  auto const jumps = find_canonical_jumps(deg, uv2v_canon, sorted2u);
  auto const e2sorted = collect_marked(jumps);
  auto const e2u = compound_maps(e2sorted, sorted2u);
  return unmap<LO>(e2u, uv2v, deg);
}

LOs find_unique(LOs const hv2v, Omega_h_Family const family, Int const high_dim,
    Int const low_dim) {
  OMEGA_H_TIME_FUNCTION;
  OMEGA_H_CHECK(high_dim > low_dim);
  OMEGA_H_CHECK(low_dim <= 2);
  OMEGA_H_CHECK(hv2v.size() % element_degree(family, high_dim, VERT) == 0);
  auto const uv2v = form_uses(hv2v, family, high_dim, low_dim);
  auto const deg = element_degree(family, low_dim, VERT);
  return find_unique_deg(deg, uv2v);
}

LOs find_unique(LOs const hv2v, Topo_type const high_type, Topo_type
    const low_type) {
  OMEGA_H_TIME_FUNCTION;
  OMEGA_H_CHECK(int(high_type) > int(low_type));
  OMEGA_H_CHECK(int(low_type) <= 3);
  OMEGA_H_CHECK(hv2v.size() % element_degree(high_type, Topo_type::vertex) == 0);
  auto const uv2v = form_uses(hv2v, high_type, low_type);
  auto const deg = element_degree(low_type, Topo_type::vertex);
  return find_unique_deg(deg, uv2v);
}

LOs form_uses(LOs const hv2v, Omega_h_Family const family, Int const high_dim,
    Int const low_dim) {
  OMEGA_H_TIME_FUNCTION;
  Int const nverts_per_high = element_degree(family, high_dim, 0);
  Int const nverts_per_low = element_degree(family, low_dim, 0);
  Int const nlows_per_high = element_degree(family, high_dim, low_dim);
  LO const nhigh = divide_no_remainder(hv2v.size(), nverts_per_high);
  LO const nuses = nhigh * nlows_per_high;
  Write<LO> uv2v(nuses * nverts_per_low);
  auto f = OMEGA_H_LAMBDA(LO h) {
    LO const h_begin = h * nverts_per_high;
    for (Int u = 0; u < nlows_per_high; ++u) {
      LO const u_begin = (h * nlows_per_high + u) * nverts_per_low;
      for (Int uv = 0; uv < nverts_per_low; ++uv) {
        uv2v[u_begin + uv] = hv2v[h_begin + element_down_template(family,
                                                high_dim, low_dim, u, uv)];
      }
    }
  };
  parallel_for(nhigh, std::move(f));
  return uv2v;
}

LOs form_uses(LOs const hv2v, Topo_type const high_type,
              Topo_type const low_type) {
  OMEGA_H_TIME_FUNCTION;
  Int const nverts_per_high = element_degree(high_type, Topo_type::vertex);
  Int const nverts_per_low = element_degree(low_type, Topo_type::vertex);
  Int const nlows_per_high = element_degree(high_type, low_type);
  LO const nhigh = divide_no_remainder(hv2v.size(), nverts_per_high);
  LO const nuses = nhigh * nlows_per_high;
  Write<LO> uv2v(nuses * nverts_per_low);
  auto f = OMEGA_H_LAMBDA(LO h) {
    LO const h_begin = h * nverts_per_high;
    for (Int u = 0; u < nlows_per_high; ++u) {
      LO const u_begin = (h * nlows_per_high + u) * nverts_per_low;
      for (Int uv = 0; uv < nverts_per_low; ++uv) {
        uv2v[u_begin + uv] = hv2v[h_begin + element_down_template(
                             int(high_type), int(low_type), u, uv)];
      }
    }
  };
  parallel_for(nhigh, std::move(f));
  return uv2v;
}

void sort_by_high_index(
    LOs const l2lh, Write<LO> const lh2h, Write<I8> const codes) {
  OMEGA_H_TIME_FUNCTION;
  LO const nl = l2lh.size() - 1;
  auto f = OMEGA_H_LAMBDA(LO const l) {
    LO const begin = l2lh[l];
    LO const end = l2lh[l + 1];
    for (LO j = begin; j < end; ++j) {
      LO k_min = j;
      GO min_h = lh2h[j];
      for (LO k = j + 1; k < end; ++k) {
        GO const h = lh2h[k];
        if (h < min_h) {
          k_min = k;
          min_h = h;
        }
      }
      swap2(lh2h[j], lh2h[k_min]);
      swap2(codes[j], codes[k_min]);
    }
  };
  parallel_for(nl, std::move(f));
}

void separate_upward_with_codes(LO const nlh, LOs const lh2hl,
    Int const nlows_per_high, Write<LO> const lh2h, Bytes const down_codes,
    Write<Byte> const codes) {
  OMEGA_H_TIME_FUNCTION;
  auto f = OMEGA_H_LAMBDA(LO lh) {
    LO const hl = lh2hl[lh];
    LO const h = hl / nlows_per_high;
    lh2h[lh] = h;
    Int const which_down = hl % nlows_per_high;
    auto const down_code = down_codes[hl];
    bool const is_flipped = code_is_flipped(down_code);
    Int const rotation = code_rotation(down_code);
    codes[lh] = make_code(is_flipped, rotation, which_down);
  };
  parallel_for(nlh, std::move(f));
}

void separate_upward_no_codes(LO const nlh, LOs const lh2hl,
    Int const nlows_per_high, Write<LO> const lh2h, Write<Byte> const codes) {
  auto f = OMEGA_H_LAMBDA(LO lh) {
    LO const hl = lh2hl[lh];
    LO const h = hl / nlows_per_high;
    lh2h[lh] = h;
    Int const which_down = hl % nlows_per_high;
    codes[lh] = make_code(false, 0, which_down);
  };
  parallel_for(nlh, std::move(f));
}

Adj invert_adj(Adj const down, Int const nlows_per_high, LO const nlows,
    Topo_type high_type, Topo_type low_type) {
  OMEGA_H_TIME_FUNCTION;
  auto const high_plural_name = dimensional_plural_name(high_type);//change code
  auto const high_singular_name = dimensional_singular_name(high_type);//change code
  auto const low_plural_name = dimensional_plural_name(low_type);
  auto const low_singular_name = dimensional_singular_name(low_type);
  auto const l2lh_name = std::string(low_plural_name) + " to " +
                         low_singular_name + " " + high_plural_name;
  auto const lh2hl_name = std::string(low_singular_name) + " " +
                          high_plural_name + " to " + high_singular_name + " " +
                          low_plural_name;
  auto const lh2h_name = std::string(low_singular_name) + " " +
                         high_plural_name + " to " + high_plural_name;
  auto const codes_name =
      std::string(low_singular_name) + " " + high_plural_name + " codes";
  auto const l2hl =
      invert_map_by_atomics(down.ab2b, nlows, l2lh_name, lh2hl_name);
  auto const l2lh = l2hl.a2ab;
  auto const lh2hl = l2hl.ab2b;
  LO const nlh = lh2hl.size();
  Read<I8> down_codes(down.codes);
  Write<LO> lh2h(nlh, lh2h_name);
  Write<I8> codes(nlh, codes_name);
  if (down_codes.exists()) {
    separate_upward_with_codes(
        nlh, lh2hl, nlows_per_high, lh2h, down_codes, codes);
  } else {
    separate_upward_no_codes(nlh, lh2hl, nlows_per_high, lh2h, codes);
  }
  sort_by_high_index(l2lh, lh2h, codes);
  return Adj(l2lh, lh2h, codes);
}

Adj invert_adj(Adj const down, Int const nlows_per_high, LO const nlows,
    Int const high_dim, Int const low_dim) {
  OMEGA_H_TIME_FUNCTION;
  auto const high_plural_name = dimensional_plural_name(high_dim);
  auto const high_singular_name = dimensional_singular_name(high_dim);
  auto const low_plural_name = dimensional_plural_name(low_dim);
  auto const low_singular_name = dimensional_singular_name(low_dim);
  auto const l2lh_name = std::string(low_plural_name) + " to " +
                         low_singular_name + " " + high_plural_name;
  auto const lh2hl_name = std::string(low_singular_name) + " " +
                          high_plural_name + " to " + high_singular_name + " " +
                          low_plural_name;
  auto const lh2h_name = std::string(low_singular_name) + " " +
                         high_plural_name + " to " + high_plural_name;
  auto const codes_name =
      std::string(low_singular_name) + " " + high_plural_name + " codes";
  auto const l2hl =
      invert_map_by_atomics(down.ab2b, nlows, l2lh_name, lh2hl_name);
  auto const l2lh = l2hl.a2ab;
  auto const lh2hl = l2hl.ab2b;
  LO const nlh = lh2hl.size();
  Read<I8> down_codes(down.codes);
  Write<LO> lh2h(nlh, lh2h_name);
  Write<I8> codes(nlh, codes_name);
  if (down_codes.exists()) {
    separate_upward_with_codes(
        nlh, lh2hl, nlows_per_high, lh2h, down_codes, codes);
  } else {
    separate_upward_no_codes(nlh, lh2hl, nlows_per_high, lh2h, codes);
  }
  sort_by_high_index(l2lh, lh2h, codes);
  return Adj(l2lh, lh2h, codes);
}

Bytes filter_parents(Parents const c2p, Int const parent_dim) {
  OMEGA_H_TIME_FUNCTION;
  Write<Byte> filter(c2p.parent_idx.size());
  auto f = OMEGA_H_LAMBDA(LO c) {
    auto const code = c2p.codes[c];
    if (amr::code_parent_dim(code) == parent_dim)
      filter[c] = 1;
    else
      filter[c] = 0;
  };
  parallel_for(c2p.parent_idx.size(), std::move(f));
  return filter;
}

Children invert_parents(
    Parents const c2p, Int const parent_dim, Int const nparent_dim_ents) {
  OMEGA_H_TIME_FUNCTION;
  auto const filter = filter_parents(c2p, parent_dim);
  auto const rc2c = collect_marked(filter);
  auto const rc2p = unmap(rc2c, c2p.parent_idx, 1);
  auto const p2rc = invert_map_by_atomics(rc2p, nparent_dim_ents);
  auto const p2pc = p2rc.a2ab;
  auto const pc2rc = p2rc.ab2b;
  auto const pc2c = unmap(pc2rc, rc2c, 1);
  auto const codes = unmap(pc2c, c2p.codes, 1);
  sort_by_high_index(p2pc, pc2c, codes);
  return Children(p2pc, pc2c, codes);
}

template <Int deg>
struct IsMatch;

// edges
template <>
struct IsMatch<2> {
  template <typename T>
  OMEGA_H_DEVICE static bool eval(Read<T> const& av2v, LO const a_begin,
      Read<T> const& bv2v, LO const b_begin, Int const which_down,
      I8* match_code) {
    if (av2v[a_begin + 1] == bv2v[b_begin + (1 - which_down)]) {
      *match_code = make_code(false, which_down, 0);
      return true;
    }
    return false;
  }
};

// triangles
template <>
struct IsMatch<3> {
  template <typename T>
  OMEGA_H_DEVICE static bool eval(Read<T> const& av2v, LO const a_begin,
      Read<T> const& bv2v, LO const b_begin, Int const which_down,
      I8* match_code) {
    if (av2v[a_begin + 1] == bv2v[b_begin + ((which_down + 1) % 3)] &&
        av2v[a_begin + 2] == bv2v[b_begin + ((which_down + 2) % 3)]) {
      *match_code = make_code(false, rotation_to_first(3, which_down), 0);
      return true;
    }
    if (av2v[a_begin + 1] == bv2v[b_begin + ((which_down + 2) % 3)] &&
        av2v[a_begin + 2] == bv2v[b_begin + ((which_down + 1) % 3)]) {
      *match_code = make_code(true, rotation_to_first(3, which_down), 0);
      return true;
    }
    return false;
  }
};

// quads
template <>
struct IsMatch<4> {
  template <typename T>
  OMEGA_H_DEVICE static bool eval(Read<T> const& av2v, LO const a_begin,
      Read<T> const& bv2v, LO const b_begin, Int const which_down,
      I8* match_code) {
    if (av2v[a_begin + 2] != bv2v[b_begin + ((which_down + 2) % 4)]) {
      return false;
    }
    if (av2v[a_begin + 1] == bv2v[b_begin + ((which_down + 1) % 4)] &&
        av2v[a_begin + 3] == bv2v[b_begin + ((which_down + 3) % 4)]) {
      *match_code = make_code(false, rotation_to_first(4, which_down), 0);
      return true;
    }
    if (av2v[a_begin + 1] == bv2v[b_begin + ((which_down + 3) % 4)] &&
        av2v[a_begin + 3] == bv2v[b_begin + ((which_down + 1) % 4)]) {
      *match_code = make_code(true, rotation_to_first(4, which_down), 0);
      return true;
    }
    return false;
  }
};

template <Int deg, typename T>
void find_matches_deg(LOs const a2fv, Read<T> const av2v,
    Read<T> const bv2v, Adj const v2b, Write<LO>* a2b_out, Write<I8>* codes_out,
    bool const allow_duplicates) {
  OMEGA_H_TIME_FUNCTION;
  LO const na = a2fv.size();
  OMEGA_H_CHECK(na * deg == av2v.size());
  LOs const v2vb = v2b.a2ab;
  LOs const vb2b = v2b.ab2b;
  Read<I8> const vb_codes = v2b.codes;
  Write<LO> a2b(na);
  Write<I8> codes(na);
  auto f = OMEGA_H_LAMBDA(LO a) {
    auto const fv = a2fv[a];
    auto const a_begin = a * deg;
    auto const vb_begin = v2vb[fv];
    auto const vb_end = v2vb[fv + 1];
    bool found = false;
    for (LO vb = vb_begin; vb < vb_end; ++vb) {
      auto const b = vb2b[vb];
      auto const vb_code = vb_codes[vb];
      auto const which_down = code_which_down(vb_code);
      auto const b_begin = b * deg;
      I8 match_code;
      if (IsMatch<deg>::eval(
              av2v, a_begin, bv2v, b_begin, which_down, &match_code)) {
        OMEGA_H_CHECK(!found); // there can't be more than one!
        a2b[a] = b;
        codes[a] = match_code;
        found = true;
        if (allow_duplicates) break;
      }
    }
    (void)found;
    OMEGA_H_CHECK(found);  // there can't be less than one!
  };
  parallel_for(na, std::move(f));
  *a2b_out = a2b;
  *codes_out = codes;
}

template <typename T>
void find_matches_ex(Int const deg, LOs const a2fv, Read<T> const av2v,
    Read<T> const bv2v, Adj const v2b, Write<LO>* a2b_out, Write<I8>* codes_out,
    bool const allow_duplicates) {
  if (deg == 2) {
    find_matches_deg<2>(
        a2fv, av2v, bv2v, v2b, a2b_out, codes_out, allow_duplicates);
  } else if (deg == 3) {
    find_matches_deg<3>(
        a2fv, av2v, bv2v, v2b, a2b_out, codes_out, allow_duplicates);
  } else if (deg == 4) {
    find_matches_deg<4>(
        a2fv, av2v, bv2v, v2b, a2b_out, codes_out, allow_duplicates);
  } else {
    Omega_h_fail("find_matches_ex called with unsupported degree %d\n", deg);
  }
}

void find_matches(Omega_h_Family const family, Int const dim, LOs const av2v,
    LOs const bv2v, Adj const v2b, Write<LO>* a2b_out, Write<I8>* codes_out) {
  OMEGA_H_CHECK(dim <= 2);
  auto const deg = element_degree(family, dim, VERT);
  auto const a2fv = get_component(av2v, deg, 0);
  find_matches_ex(deg, a2fv, av2v, bv2v, v2b, a2b_out, codes_out);
}

void find_matches(Topo_type const ent_type, LOs const av2v,
    LOs const bv2v, Adj const v2b, Write<LO>* a2b_out, Write<I8>* codes_out) {
  OMEGA_H_CHECK(int(ent_type) <= 3);
  auto const deg = element_degree(ent_type, Topo_type::vertex);
  auto const a2fv = get_component(av2v, deg, 0);
  find_matches_ex(deg, a2fv, av2v, bv2v, v2b, a2b_out, codes_out);
}

Adj reflect_down(LOs const hv2v, LOs const lv2v, Adj const v2l,
    Omega_h_Family const family, Int const high_dim, Int const low_dim) {
  ScopedTimer timer("reflect_down(v2l)");
  LOs const uv2v = form_uses(hv2v, family, high_dim, low_dim);
  Write<LO> hl2l;
  Write<I8> codes;
  find_matches(family, low_dim, uv2v, lv2v, v2l, &hl2l, &codes);
  return Adj(read(hl2l), read(codes));
}

Adj reflect_down(LOs const hv2v, LOs const lv2v, Adj const v2l,
    Topo_type const high_type, Topo_type const low_type) {
  ScopedTimer timer("reflect_down_mix(v2l)");
  LOs const uv2v = form_uses(hv2v, high_type, low_type);
  Write<LO> hl2l;
  Write<I8> codes;
  find_matches(low_type, uv2v, lv2v, v2l, &hl2l, &codes);
  return Adj(read(hl2l), read(codes));
}

Adj reflect_down(LOs const hv2v, LOs const lv2v, Omega_h_Family const family,
    LO const nv, Int const high_dim, Int const low_dim) {
  ScopedTimer timer("reflect_down(nv)");
  auto const nverts_per_low = element_degree(family, low_dim, 0);
  auto const l2v = Adj(lv2v);
  auto const v2l = invert_adj(l2v, nverts_per_low, nv, high_dim, low_dim);
  return reflect_down(hv2v, lv2v, v2l, family, high_dim, low_dim);
}

Adj transit(Adj const h2m, Adj const m2l, Omega_h_Family const family,
    Int const high_dim, Int const low_dim) {
  OMEGA_H_TIME_FUNCTION;
  auto const high_singular_name = dimensional_singular_name(high_dim);
  auto const low_plural_name = dimensional_plural_name(low_dim);
  auto const hl2l_name = std::string(high_singular_name) + " " +
                         low_plural_name + " to " + low_plural_name;
  OMEGA_H_CHECK(3 >= high_dim);
  auto const mid_dim = low_dim + 1;
  OMEGA_H_CHECK(high_dim > mid_dim);
  OMEGA_H_CHECK(low_dim == 1 || low_dim == 0);
  auto const hm2m = h2m.ab2b;
  auto const m2hm_codes = h2m.codes;
  auto const ml2l = m2l.ab2b;
  auto const ml_codes = m2l.codes;
  auto const nmids_per_high = element_degree(family, high_dim, mid_dim);
  auto const nlows_per_mid = element_degree(family, mid_dim, low_dim);
  auto const nlows_per_high = element_degree(family, high_dim, low_dim);
  auto const nhighs = divide_no_remainder(hm2m.size(), nmids_per_high);
  Write<LO> hl2l(nhighs * nlows_per_high, hl2l_name);
  Write<I8> codes;
  /* codes only need to be created when transiting region->face + face->edge =
     region->edge. any other transit has vertices as its destination, and
     vertices have no orientation/alignment */
  if (low_dim == 1) {
    auto const codes_name =
        std::string(high_singular_name) + " " + low_plural_name + " codes";
    codes = Write<I8>(hl2l.size(), codes_name);
  }
  auto f = OMEGA_H_LAMBDA(LO h) {
    auto const hl_begin = h * nlows_per_high;
    auto const hm_begin = h * nmids_per_high;
    for (Int hl = 0; hl < nlows_per_high; ++hl) {
      auto const ut = element_up_template(family, high_dim, low_dim, hl, 0);
      auto const hm = ut.up;
      auto const hml = ut.which_down;
      auto const m = hm2m[hm_begin + hm];
      auto const m2hm_code = m2hm_codes[hm_begin + hm];
      auto const hm2m_code = invert_alignment(nlows_per_mid, m2hm_code);
      auto const ml = align_index(nlows_per_mid, low_dim, hml, hm2m_code);
      auto const ml_begin = m * nlows_per_mid;
      auto const l = ml2l[ml_begin + ml];
      // safety check for duplicates.
      // remove after this code is heavily exercised (or don't)
      for (Int hhl2 = 0; hhl2 < hl; ++hhl2) {
        OMEGA_H_CHECK(l != hl2l[hl_begin + hhl2]);
      }
      hl2l[hl_begin + hl] = l;
      if (low_dim == 1) {
        /* all we are determining here is whether the edge is pointed
           in or against the direction of the "canonical edge" as
           defined by the element's template.
           this is a bitwise XOR of several flips along the way */
        auto const region_face_code = hm2m_code;
        auto const face_edge_code = ml_codes[ml_begin + ml];
        auto const region_face_flipped = code_is_flipped(region_face_code);
        auto const face_edge_flipped = bool(code_rotation(face_edge_code) == 1);
        auto const canon_flipped = ut.is_flipped;
        bool const region_edge_flipped =
            region_face_flipped ^ face_edge_flipped ^ canon_flipped;
        codes[hl_begin + hl] = make_code(false, region_edge_flipped, 0);
      }
    }
  };
  parallel_for(nhighs, std::move(f));
  if (low_dim == 1) return Adj(hl2l, codes);
  return Adj(hl2l);
}

Graph verts_across_edges(Adj const e2v, Adj const v2e) {
  OMEGA_H_TIME_FUNCTION;
  auto const ev2v = e2v.ab2b;
  auto const v2ve = v2e.a2ab;
  auto const ve2e = v2e.ab2b;
  auto const v2ve_codes = v2e.codes;
  auto const v2vv = v2ve;
  Write<LO> vv2v(ve2e.size());
  auto f = OMEGA_H_LAMBDA(LO ve) {
    auto const vv = ve;
    auto const e = ve2e[ve];
    auto const v2ve_code = v2ve_codes[ve];
    auto const eev = code_which_down(v2ve_code);
    auto const v = ev2v[e * 2 + (1 - eev)];
    vv2v[vv] = v;
  };
  parallel_for(vv2v.size(), std::move(f));
  return Adj(v2vv, vv2v);
}

Adj transit(Adj const h2m, Adj const m2l, 
    Topo_type const high_type, Topo_type const low_type, Topo_type const mid_type) {
  OMEGA_H_TIME_FUNCTION;
  auto const high_singular_name = dimensional_singular_name(high_type);
  auto const low_plural_name = dimensional_plural_name(low_type);
  auto const hl2l_name = std::string(high_singular_name) + " " +
                         low_plural_name + " to " + low_plural_name;
  OMEGA_H_CHECK(7 >= int(high_type));
  OMEGA_H_CHECK(int(high_type) > int(mid_type));
  OMEGA_H_CHECK(int(low_type) == 1 || int(low_type) == 0);
  auto const hm2m = h2m.ab2b;
  auto const m2hm_codes = h2m.codes;
  auto const ml2l = m2l.ab2b;
  auto const ml_codes = m2l.codes;
  auto const nmids_per_high = element_degree(high_type, mid_type);
  auto const nlows_per_mid = element_degree(mid_type, low_type);
  auto const nlows_per_high = element_degree(high_type, low_type);
  auto const nhighs = divide_no_remainder(hm2m.size(), nmids_per_high);
  Write<LO> hl2l(nhighs * nlows_per_high, hl2l_name);
  Write<I8> codes;
  if (int(low_type) == 1) {
    auto const codes_name =
        std::string(high_singular_name) + " " + low_plural_name + " codes";
    codes = Write<I8>(hl2l.size(), codes_name);
  }
  auto f = OMEGA_H_LAMBDA(LO h) {
    auto const hl_begin = h * nlows_per_high;
    auto const hm_begin = h * nmids_per_high;
    for (Int hl = 0; hl < nlows_per_high; ++hl) {
      auto const ut = element_up_template(int(high_type), int(low_type), hl, 0);//add wedge and pyramid templates
      auto const hm = ut.up;
      auto const hml = ut.which_down;
      auto const m = hm2m[hm_begin + hm];
      auto const m2hm_code = m2hm_codes[hm_begin + hm];
      auto const hm2m_code = invert_alignment(nlows_per_mid, m2hm_code);
      auto const ml = align_index(nlows_per_mid, int(low_type), hml, hm2m_code);
      auto const ml_begin = m * nlows_per_mid;
      auto const l = ml2l[ml_begin + ml];
      for (Int hhl2 = 0; hhl2 < hl; ++hhl2) {
        OMEGA_H_CHECK(l != hl2l[hl_begin + hhl2]);
      }
      hl2l[hl_begin + hl] = l;
      if (int(low_type) == 1) {
        auto const region_face_code = hm2m_code;
        auto const face_edge_code = ml_codes[ml_begin + ml];
        auto const region_face_flipped = code_is_flipped(region_face_code);
        auto const face_edge_flipped = bool(code_rotation(face_edge_code) == 1);
        auto const canon_flipped = ut.is_flipped;
        bool const region_edge_flipped =
            region_face_flipped ^ face_edge_flipped ^ canon_flipped;
        codes[hl_begin + hl] = make_code(false, region_edge_flipped, 0);
      }
    }
  };
  parallel_for(nhighs, std::move(f));
  if (int(low_type) == 1) return Adj(hl2l, codes);
  return Adj(hl2l);
}

Graph edges_across_tris(Adj const f2e, Adj const e2f) {
  OMEGA_H_TIME_FUNCTION;
  auto const fe2e = f2e.ab2b;
  auto const e2ef = e2f.a2ab;
  auto const ef2f = e2f.ab2b;
  auto const e2ef_codes = e2f.codes;
  auto const ne = e2ef.size() - 1;
  auto const e2ef_degrees = get_degrees(e2ef);
  auto const e2ee_degrees = multiply_each_by(e2ef_degrees, 2);
  auto const e2ee = offset_scan(e2ee_degrees, "edges to edge edges");
  auto const nee = e2ee.last();
  Write<LO> ee2e(nee, "edge edges to edges");
  auto functor = OMEGA_H_LAMBDA(LO e) {
    auto const ef_begin = e2ef[e];
    auto const ef_end = e2ef[e + 1];
    auto const neef = ef_end - ef_begin;
    auto const ee_begin = e2ee[e];
    for (Int eef = 0; eef < neef; ++eef) {
      auto const ef = ef_begin + eef;
      auto const f = ef2f[ef];
      auto const e2ef_code = e2ef_codes[ef];
      auto const ffe = code_which_down(e2ef_code);
      auto const e1 = fe2e[f * 3 + ((ffe + 1) % 3)];
      auto const e2 = fe2e[f * 3 + ((ffe + 2) % 3)];
      ee2e[ee_begin + eef * 2 + 0] = e1;
      ee2e[ee_begin + eef * 2 + 1] = e2;
    }
  };
  parallel_for(ne, std::move(functor));
  return Adj(e2ee, ee2e);
}

Graph edges_across_tets(Adj const r2e, Adj const e2r) {
  OMEGA_H_TIME_FUNCTION;
  auto const re2e = r2e.ab2b;
  auto const e2er = e2r.a2ab;
  auto const er2r = e2r.ab2b;
  auto const e2er_codes = e2r.codes;
  auto const ne = e2er.size() - 1;
  auto const e2ee = e2er;
  Write<LO> ee2e(er2r.size(), "edge edges to edges");
  auto f = OMEGA_H_LAMBDA(LO e) {
    auto const er_begin = e2er[e];
    auto const er_end = e2er[e + 1];
    for (auto er = er_begin; er < er_end; ++er) {
      auto const r = er2r[er];
      auto const e2er_code = e2er_codes[er];
      auto const rre = code_which_down(e2er_code);
      auto const rre_opp = simplex_opposite_template(REGION, EDGE, rre);
      auto const re_begin = r * 6;
      auto const e_opp = re2e[re_begin + rre_opp];
      auto const ee = er;
      ee2e[ee] = e_opp;
    }
  };
  parallel_for(ne, std::move(f));
  return Adj(e2ee, ee2e);
}

Graph elements_across_sides(Int const dim, Adj const elems2sides,
    Adj const sides2elems, Read<I8> const side_is_exposed) {
  OMEGA_H_TIME_FUNCTION;
  auto const elem_side2side = elems2sides.ab2b;
  auto const side2side_elems = sides2elems.a2ab;
  auto const side_elem2elem = sides2elems.ab2b;
  Int const nsides_per_elem = dim + 1;
  auto const nelems =
      divide_no_remainder(elem_side2side.size(), nsides_per_elem);
  Write<LO> degrees(nelems);
  auto count = OMEGA_H_LAMBDA(LO elem) {
    auto const begin = elem * nsides_per_elem;
    auto const end = begin + nsides_per_elem;
    Int n = 0;
    for (auto elem_side = begin; elem_side < end; ++elem_side) {
      auto const side = elem_side2side[elem_side];
      if (!side_is_exposed[side]) {
        OMEGA_H_CHECK(side2side_elems[side + 1] - side2side_elems[side] == 2);
        ++n;
      }
    }
    degrees[elem] = n;
  };
  parallel_for(nelems, std::move(count));
  auto const elem2elem_elems = offset_scan(LOs(degrees));
  auto const nelem_elems = elem2elem_elems.last();
  Write<LO> elem_elem2elem(nelem_elems);
  auto fill = OMEGA_H_LAMBDA(LO elem) {
    auto const begin = elem * nsides_per_elem;
    auto const end = begin + nsides_per_elem;
    LO elem_elem = elem2elem_elems[elem];
    for (auto elem_side = begin; elem_side < end; ++elem_side) {
      auto const side = elem_side2side[elem_side];
      if (!side_is_exposed[side]) {
        auto side_elem = side2side_elems[side];
        if (side_elem2elem[side_elem] == elem) ++side_elem;
        elem_elem2elem[elem_elem] = side_elem2elem[side_elem];
        ++elem_elem;
      }
    }
  };
  parallel_for(nelems, std::move(fill));
  return Graph(elem2elem_elems, elem_elem2elem);
}

#define INST(T)                                                                \
  template Read<I8> get_codes_to_canonical(Int deg, Read<T> ev2v);             \
  template void find_matches_ex(Int deg, LOs a2fv, Read<T> av2v, Read<T> bv2v, \
      Adj v2b, Write<LO>* a2b_out, Write<I8>* codes_out, bool);
INST(LO)
INST(GO)
#undef INST

}  // namespace Omega_h
