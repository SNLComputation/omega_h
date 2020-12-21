#include "Omega_h_class.hpp"

#include "Omega_h_array_ops.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mark.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_surface.hpp"

namespace Omega_h {

void classify_sides_by_exposure(Mesh* mesh, Read<I8> side_is_exposed) {
  auto dim = mesh->dim();
  auto ns = mesh->nents(dim - 1);
  Write<I8> class_dim(ns);
  auto f = OMEGA_H_LAMBDA(LO s) {
    class_dim[s] = static_cast<I8>(dim - side_is_exposed[s]);
  };
  parallel_for(ns, f, "classify_sides_by_exposure");
  mesh->add_tag<I8>(dim - 1, "class_dim", 1, class_dim);
}

void classify_hinges_by_sharpness(
    Mesh* mesh, Read<I8> hinge_is_exposed, Read<I8> hinge_is_sharp) {
  auto dim = mesh->dim();
  auto nh = mesh->nents(dim - 2);
  Write<I8> class_dim(nh);
  auto f = OMEGA_H_LAMBDA(LO h) {
    class_dim[h] =
        static_cast<I8>(dim - hinge_is_exposed[h] - hinge_is_sharp[h]);
  };
  parallel_for(nh, f, "classify_hinges_by_sharpness");
  mesh->add_tag<I8>(dim - 2, "class_dim", 1, class_dim);
}

void classify_elements(Mesh* mesh) {
  mesh->add_tag<I8>(mesh->dim(), "class_dim", 1,

      Read<I8>(mesh->nelems(), static_cast<I8>(mesh->dim())));
}

void classify_by_angles(Mesh* mesh, Real sharp_angle) {
  OMEGA_H_CHECK(mesh->family() == OMEGA_H_SIMPLEX);
  auto dim = mesh->dim();
  classify_elements(mesh);
  auto side_is_exposed = mark_exposed_sides(mesh);
  classify_sides_by_exposure(mesh, side_is_exposed);
  if (dim == 1) return;
  auto hinge_is_exposed = mark_down(mesh, dim - 1, dim - 2, side_is_exposed);
  auto surf_side2side = collect_marked(side_is_exposed);
  auto surf_side_normals = get_side_vectors(mesh, surf_side2side);
  auto surf_hinge2hinge = collect_marked(hinge_is_exposed);
  auto nsurf_hinges = surf_hinge2hinge.size();
  auto nsides = mesh->nents(dim - 1);
  auto side2surf_side = invert_injective_map(surf_side2side, nsides);
  auto surf_hinge_angles = get_hinge_angles(
      mesh, surf_side_normals, surf_hinge2hinge, side2surf_side);
  auto nhinges = mesh->nents(dim - 2);
  Write<I8> hinge_is_sharp(nhinges, 0);
  auto f = OMEGA_H_LAMBDA(LO surf_hinge) {
    LO hinge = surf_hinge2hinge[surf_hinge];
    hinge_is_sharp[hinge] = (surf_hinge_angles[surf_hinge] >= sharp_angle);
  };
  parallel_for(nsurf_hinges, f, "classify_by_angles(is_sharp)");
  classify_hinges_by_sharpness(mesh, hinge_is_exposed, hinge_is_sharp);
  if (dim == 2) return;
  finalize_classification(mesh);
}

static bool has_any_ids(Mesh* mesh) {
  for (Int dim = 0; dim <= mesh->dim(); ++dim) {
    if (mesh->has_tag(dim, "class_id")) return true;
  }
  return false;
}

static void remove_all_ids(Mesh* mesh) {
  for (Int dim = 0; dim <= mesh->dim(); ++dim) {
    mesh->remove_tag(dim, "class_id");
  }
}

template <typename T>
static Write<T> copy_and_remove_or_default(
    Mesh* mesh, Int dim, std::string const& name, T def_val) {
  if (mesh->has_tag(dim, name)) {
    auto a = deep_copy(mesh->get_array<T>(dim, name));
    mesh->remove_tag(dim, name);
    return a;
  } else {
    return Write<T>(mesh->nents(dim), def_val);
  }
}

void project_classification(Mesh* mesh, Int ent_dim, Write<I8> class_dim,
    Write<ClassId> class_id, bool relaxed) {
  auto l2h = mesh->ask_up(ent_dim, ent_dim + 1);
  auto l2lh = l2h.a2ab;
  auto lh2h = l2h.ab2b;
  auto high_class_dim = mesh->get_array<I8>(ent_dim + 1, "class_dim");
  auto high_class_id = mesh->get_array<ClassId>(ent_dim + 1, "class_id");
  auto f = OMEGA_H_LAMBDA(LO l) {
    Int best_dim = class_dim[l];
    auto best_id = class_id[l];
    Int nadj = 0;
    for (auto lh = l2lh[l]; lh < l2lh[l + 1]; ++lh) {
      auto h = lh2h[lh];
      Int high_dim = high_class_dim[h];
      auto high_id = high_class_id[h];
      if (high_dim < best_dim) {
        best_dim = high_dim;
        best_id = high_id;
        nadj = 1;
      } else if (high_dim == best_dim) {
        if (high_id != best_id && (!relaxed || high_id != -1)) {
          if (best_id == -1) {
            best_id = high_id;
            ++nadj;
          } else {
            --best_dim;
            best_id = -1;
            nadj = 0;
          }
        } else {
          ++nadj;
        }
      }
    }
    if ((nadj != 2 && best_dim == ent_dim + 1) ||
        (nadj < 2 && best_dim > ent_dim + 1)) {
      best_dim = ent_dim;
      best_id = -1;
    }
    class_dim[l] = static_cast<I8>(best_dim);
    class_id[l] = best_id;
  };
  parallel_for(mesh->nents(ent_dim), f, "project_classification");
}

void project_classification(Mesh* mesh, Int ent_dim, bool relaxed) {
  OMEGA_H_CHECK(mesh->comm()->size() == 1);
  auto class_ids_w = deep_copy(mesh->get_array<ClassId>(ent_dim, "class_id"));
  auto class_dims_w = deep_copy(mesh->get_array<I8>(ent_dim, "class_dim"));
  project_classification(mesh, ent_dim, class_dims_w, class_ids_w, relaxed);
  mesh->set_tag(ent_dim, "class_id", Read<ClassId>(class_ids_w));
  mesh->set_tag(ent_dim, "class_dim", Read<I8>(class_dims_w));
}

void finalize_classification(Mesh* mesh) {
  bool had_ids = has_any_ids(mesh);
  for (Int ent_dim = mesh->dim(); ent_dim >= VERT; --ent_dim) {
    OMEGA_H_CHECK(mesh->owners_have_all_upward(ent_dim));
    auto class_dim_w = copy_and_remove_or_default<I8>(
        mesh, ent_dim, "class_dim", I8(mesh->dim()));
    auto class_id_w =
        copy_and_remove_or_default<ClassId>(mesh, ent_dim, "class_id", -1);
    if (ent_dim < mesh->dim()) {
      project_classification(mesh, ent_dim, class_dim_w, class_id_w);
    }
    auto class_dim = mesh->sync_array(ent_dim, read(class_dim_w), 1);
    auto class_id = mesh->sync_array(ent_dim, read(class_id_w), 1);
    mesh->add_tag<I8>(ent_dim, "class_dim", 1, class_dim);
    mesh->add_tag<ClassId>(ent_dim, "class_id", 1, class_id);
  }
  if (!had_ids) remove_all_ids(mesh);
}

void classify_equal_order(
    Mesh* mesh, Int ent_dim, LOs eqv2v, Read<ClassId> eq_class_ids) {
  LOs eq2e;
  if (ent_dim == mesh->dim()) {
    /* assuming elements were constructed in the same order ! */
    eq2e = LOs(mesh->nelems(), 0, 1);
  } else if (ent_dim == VERT) {
    eq2e = eqv2v;
  } else {
    Write<LO> eq2e_w;
    Write<I8> codes;
    auto ev2v = mesh->ask_verts_of(ent_dim);
    auto v2e = mesh->ask_up(VERT, ent_dim);
    find_matches(mesh->family(), ent_dim, eqv2v, ev2v, v2e, &eq2e_w, &codes);
    eq2e = eq2e_w;
  }
  auto neq = eqv2v.size() / (ent_dim + 1);
  auto eq_class_dim = Read<I8>(neq, I8(ent_dim));
  auto class_dim =
      map_onto(eq_class_dim, eq2e, mesh->nents(ent_dim), I8(mesh->dim()), 1);
  auto class_id = map_onto(eq_class_ids, eq2e, mesh->nents(ent_dim), -1, 1);
  mesh->add_tag<I8>(ent_dim, "class_dim", 1, class_dim);
  mesh->add_tag<ClassId>(ent_dim, "class_id", 1, class_id);
}

/* TODO this function is included as work in progress and is not tested */
void classify_equal_order(
    Mesh* mesh, Topo_type ent_type, LOs eqv2v, Read<ClassId> eq_class_ids) {
  LOs eq2e;
  if (ent_type == Topo_type::vertex) {
    eq2e = eqv2v;
  } 
  else {
    Write<LO> eq2e_w;
    Write<I8> codes;
    auto ev2v = mesh->ask_verts_of(ent_type);
    auto v2e = mesh->ask_up(Topo_type::vertex, ent_type);
    find_matches(ent_type, eqv2v, ev2v, v2e, &eq2e_w, &codes);
    eq2e = eq2e_w;
  }

  auto ent_dim = mesh->ent_dim(ent_type);

  auto neq = eqv2v.size() / (ent_dim + 1);
  auto eq_class_dim = Read<I8>(neq, I8(ent_dim));
  auto class_dim =
      map_onto(eq_class_dim, eq2e, mesh->nents(ent_type), I8(mesh->dim()), 1);
  auto class_id = map_onto(eq_class_ids, eq2e, mesh->nents(ent_type), -1, 1);
  mesh->add_tag<I8>(ent_type, "class_dim", 1, class_dim);
  mesh->add_tag<ClassId>(ent_type, "class_id", 1, class_id);
}
/* */

}  // end namespace Omega_h
