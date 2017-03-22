#include "Omega_h_metric.hpp"

#include <iostream>

#include "Omega_h_array_ops.hpp"
#include "host_few.hpp"
#include "project.hpp"
#include "Omega_h_simplex.hpp"
#include "size.hpp"

namespace Omega_h {

Int get_metric_dim(Int ncomps) {
  for (Int i = 1; i <= 3; ++i)
    if (ncomps == symm_ncomps(i)) return i;
  NORETURN(Int());
}

Int get_metrics_dim(LO nmetrics, Reals metrics) {
  CHECK(metrics.size() % nmetrics == 0);
  auto ncomps = metrics.size() / nmetrics;
  return get_metric_dim(ncomps);
}

Int get_metric_dim(Mesh* mesh) {
  auto ncomps = mesh->get_tagbase(VERT, "metric")->ncomps();
  return get_metric_dim(ncomps);
}

template <Int dim>
static Reals clamp_metrics_dim(
    LO nmetrics, Reals metrics, Real h_min, Real h_max) {
  auto out = Write<Real>(nmetrics * symm_ncomps(dim));
  auto f = LAMBDA(LO i) {
    auto m = get_symm<dim>(metrics, i);
    m = clamp_metric(m, h_min, h_max);
    set_symm(out, i, m);
  };
  parallel_for(nmetrics, f);
  return out;
}

Reals clamp_metrics(LO nmetrics, Reals metrics, Real h_min, Real h_max) {
  auto dim = get_metrics_dim(nmetrics, metrics);
  if (dim == 3) return clamp_metrics_dim<3>(nmetrics, metrics, h_min, h_max);
  if (dim == 2) return clamp_metrics_dim<2>(nmetrics, metrics, h_min, h_max);
  if (dim == 1) return clamp_metrics_dim<1>(nmetrics, metrics, h_min, h_max);
  NORETURN(Reals());
}

template <Int mdim, Int edim>
static Reals mident_metrics_tmpl(Mesh* mesh, LOs a2e, Reals v2m) {
  auto na = a2e.size();
  Write<Real> out(na * symm_ncomps(mdim));
  auto ev2v = mesh->ask_verts_of(edim);
  auto f = LAMBDA(LO a) {
    auto e = a2e[a];
    auto v = gather_verts<edim + 1>(ev2v, e);
    auto ms = gather_symms<edim + 1, mdim>(v2m, v);
    auto m = average_metric(ms);
    set_symm(out, a, m);
  };
  parallel_for(na, f);
  return out;
}

Reals get_mident_metrics(Mesh* mesh, Int ent_dim, LOs entities, Reals v2m) {
  auto metrics_dim = get_metrics_dim(mesh->nverts(), v2m);
  if (metrics_dim == 3 && ent_dim == 3) {
    return mident_metrics_tmpl<3, 3>(mesh, entities, v2m);
  }
  if (metrics_dim == 3 && ent_dim == 1) {
    return mident_metrics_tmpl<3, 1>(mesh, entities, v2m);
  }
  if (metrics_dim == 2 && ent_dim == 2) {
    return mident_metrics_tmpl<2, 2>(mesh, entities, v2m);
  }
  if (metrics_dim == 2 && ent_dim == 1) {
    return mident_metrics_tmpl<2, 1>(mesh, entities, v2m);
  }
  if (metrics_dim == 1 && ent_dim == 3) {
    return mident_metrics_tmpl<1, 3>(mesh, entities, v2m);
  }
  if (metrics_dim == 1 && ent_dim == 2) {
    return mident_metrics_tmpl<1, 2>(mesh, entities, v2m);
  }
  if (metrics_dim == 1 && ent_dim == 1) {
    return mident_metrics_tmpl<1, 1>(mesh, entities, v2m);
  }
  NORETURN(Reals());
}

Reals interpolate_between_metrics(LO nmetrics, Reals a, Reals b, Real t) {
  auto log_a = linearize_metrics(nmetrics, a);
  auto log_b = linearize_metrics(nmetrics, b);
  auto log_c = interpolate_between(log_a, log_b, t);
  return delinearize_metrics(nmetrics, log_c);
}

template <Int dim>
Reals linearize_metrics_dim(Reals metrics) {
  auto n = metrics.size() / symm_ncomps(dim);
  auto out = Write<Real>(n * symm_ncomps(dim));
  auto f = LAMBDA(LO i) {
    set_symm(out, i, linearize_metric(get_symm<dim>(metrics, i)));
  };
  parallel_for(n, f);
  return out;
}

template <Int dim>
Reals delinearize_metrics_dim(Reals lms) {
  auto n = lms.size() / symm_ncomps(dim);
  auto out = Write<Real>(n * symm_ncomps(dim));
  auto f = LAMBDA(LO i) {
    set_symm(out, i, delinearize_metric(get_symm<dim>(lms, i)));
  };
  parallel_for(n, f);
  return out;
}

Reals linearize_metrics(LO nmetrics, Reals metrics) {
  auto dim = get_metrics_dim(nmetrics, metrics);
  if (dim == 3) return linearize_metrics_dim<3>(metrics);
  if (dim == 2) return linearize_metrics_dim<2>(metrics);
  if (dim == 1) return linearize_metrics_dim<1>(metrics);
  NORETURN(Reals());
}

Reals delinearize_metrics(LO nmetrics, Reals linear_metrics) {
  auto dim = get_metrics_dim(nmetrics, linear_metrics);
  if (dim == 3) return delinearize_metrics_dim<3>(linear_metrics);
  if (dim == 2) return delinearize_metrics_dim<2>(linear_metrics);
  if (dim == 1) return delinearize_metrics_dim<1>(linear_metrics);
  NORETURN(Reals());
}

template <Int dim>
static HostFew<Reals, dim> axes_from_metrics_dim(Reals metrics) {
  CHECK(metrics.size() % symm_ncomps(dim) == 0);
  auto n = metrics.size() / symm_ncomps(dim);
  HostFew<Write<Real>, dim> w;
  for (Int i = 0; i < dim; ++i) w[i] = Write<Real>(n * dim);
  auto f = LAMBDA(LO i) {
    auto md = decompose_metric(get_symm<dim>(metrics, i));
    for (Int j = 0; j < dim; ++j) set_vector(w[j], i, md.q[j] * md.l[j]);
  };
  parallel_for(n, f);
  HostFew<Reals, dim> r;
  for (Int i = 0; i < dim; ++i) r[i] = Reals(w[i]);
  return r;
}

template <Int dim>
static void axes_from_metric_field_dim(Mesh* mesh,
    std::string const& metric_name, std::string const& output_prefix) {
  auto metrics = mesh->get_array<Real>(VERT, metric_name);
  auto axes = axes_from_metrics_dim<dim>(metrics);
  for (Int i = 0; i < dim; ++i) {
    mesh->add_tag(VERT, output_prefix + '_' + std::to_string(i), dim, axes[i]);
  }
}

void axes_from_metric_field(Mesh* mesh, std::string const& metric_name,
    std::string const& axis_prefix) {
  if (mesh->dim() == 3) {
    axes_from_metric_field_dim<3>(mesh, metric_name, axis_prefix);
    return;
  }
  if (mesh->dim() == 2) {
    axes_from_metric_field_dim<2>(mesh, metric_name, axis_prefix);
    return;
  }
  NORETURN();
}

/* A Hessian-based anisotropic size field, from
 * Alauzet's tech report:
 *
 * F. Alauzet, P.J. Frey, Estimateur d'erreur geometrique
 * et metriques anisotropes pour l'adaptation de maillage.
 * Partie I: aspects theoriques,
 * RR-4759, INRIA Rocquencourt, 2003.
 */

template <Int dim>
static INLINE Matrix<dim, dim> metric_from_hessian(
    Matrix<dim, dim> hessian, Real eps) {
  auto ed = decompose_eigen(hessian);
  auto r = ed.q;
  auto l = ed.l;
  constexpr auto c_num = square(dim);
  constexpr auto c_denom = 2 * square(dim + 1);
  decltype(l) tilde_l;
  for (Int i = 0; i < dim; ++i) {
    tilde_l[i] = (c_num * fabs(l[i])) / (c_denom * eps);
  }
  return compose_eigen(r, tilde_l);
}

template <Int dim>
static Reals metric_from_hessians_dim(Reals hessians, Real eps) {
  auto ncomps = symm_ncomps(dim);
  CHECK(hessians.size() % ncomps == 0);
  auto n = hessians.size() / ncomps;
  auto out = Write<Real>(n * ncomps);
  auto f = LAMBDA(LO i) {
    auto hess = get_symm<dim>(hessians, i);
    auto m = metric_from_hessian(hess, eps);
    set_symm(out, i, m);
  };
  parallel_for(n, f);
  return out;
}

Reals metric_from_hessians(Int dim, Reals hessians, Real eps) {
  CHECK(eps > 0.0);
  if (dim == 3) return metric_from_hessians_dim<3>(hessians, eps);
  if (dim == 2) return metric_from_hessians_dim<2>(hessians, eps);
  NORETURN(Reals());
}

/* gradation limiting code: */

template <Int mesh_dim, Int metric_dim>
static Reals limit_gradation_once_tmpl(
    Mesh* mesh, Reals values, Real max_rate) {
  auto v2v = mesh->ask_star(VERT);
  auto coords = mesh->coords();
  auto out = Write<Real>(mesh->nverts() * symm_ncomps(metric_dim));
  auto f = LAMBDA(LO v) {
    auto m = get_symm<metric_dim>(values, v);
    auto x = get_vector<mesh_dim>(coords, v);
    for (auto vv = v2v.a2ab[v]; vv < v2v.a2ab[v + 1]; ++vv) {
      auto av = v2v.ab2b[vv];
      auto am = get_symm<metric_dim>(values, av);
      auto ax = get_vector<mesh_dim>(coords, av);
      auto vec = ax - x;
      auto metric_dist = metric_length(m, vec);
      auto decomp = decompose_metric(m);
      decomp.l = decomp.l * (1.0 + metric_dist * max_rate);
      auto limiter = compose_metric(decomp.q, decomp.l);
      auto limited = intersect_metrics(m, limiter);
      m = limited;
    }
    set_symm(out, v, m);
  };
  parallel_for(mesh->nverts(), f);
  values = Reals(out);
  values = mesh->sync_array(VERT, values, symm_ncomps(metric_dim));
  return values;
}

static Reals limit_gradation_once(Mesh* mesh, Reals values, Real max_rate) {
  auto metric_dim = get_metrics_dim(mesh->nverts(), values);
  if (mesh->dim() == 3 && metric_dim == 3) {
    return limit_gradation_once_tmpl<3, 3>(mesh, values, max_rate);
  }
  if (mesh->dim() == 2 && metric_dim == 2) {
    return limit_gradation_once_tmpl<2, 2>(mesh, values, max_rate);
  }
  if (mesh->dim() == 3 && metric_dim == 1) {
    return limit_gradation_once_tmpl<3, 1>(mesh, values, max_rate);
  }
  if (mesh->dim() == 2 && metric_dim == 1) {
    return limit_gradation_once_tmpl<2, 1>(mesh, values, max_rate);
  }
  NORETURN(Reals());
}

Reals limit_metric_gradation(
    Mesh* mesh, Reals values, Real max_rate, Real tol) {
  CHECK(mesh->owners_have_all_upward(VERT));
  CHECK(max_rate > 0.0);
  auto comm = mesh->comm();
  Reals values2 = values;
  Int i = 0;
  do {
    values = values2;
    values2 = limit_gradation_once(mesh, values, max_rate);
    ++i;
    if (can_print(mesh) && i > 40) {
      std::cout << "warning: gradation limiting is up to step " << i << '\n';
    }
  } while (!comm->reduce_and(are_close(values, values2, tol)));
  if (can_print(mesh)) {
    std::cout << "limited gradation in " << i << " steps\n";
  }
  return values2;
}

Reals project_metrics(Mesh* mesh, Reals e2m) {
  auto e_linear = linearize_metrics(mesh->nelems(), e2m);
  auto v_linear = project_by_average(mesh, e_linear);
  return delinearize_metrics(mesh->nverts(), v_linear);
}

Reals smooth_metric_once(Mesh* mesh, Reals v2m) {
  auto e2e = LOs(mesh->nelems(), 0, 1);
  return project_metrics(mesh, get_mident_metrics(mesh, mesh->dim(), e2e, v2m));
}

/* Micheletti, S., and S. Perotto.
 * "Anisotropic adaptation via a Zienkiewicz–Zhu error estimator
 *  for 2D elliptic problems."
 * Numerical Mathematics and Advanced Applications 2009.
 * Springer Berlin Heidelberg, 2010. 645-653.
 *
 * Farrell, P. E., S. Micheletti, and S. Perotto.
 * "An anisotropic Zienkiewicz–Zhu‐type error estimator for 3D applications."
 * International journal for numerical methods in engineering
 * 85.6 (2011): 671-692.
 */

template <Int dim>
Reals get_aniso_zz_metric_dim(
    Mesh* mesh, Reals elem_gradients, Real error_bound, Real max_size) {
  CHECK(mesh->have_all_upward());
  constexpr auto nverts_per_elem = dim + 1;
  auto elem_verts2vert = mesh->ask_elem_verts();
  auto verts2elems = mesh->ask_up(VERT, dim);
  constexpr auto max_elems_per_patch =
      AvgDegree<dim, 0, dim>::value * nverts_per_elem;
  auto elems2volume = measure_elements_real(mesh);
  auto nglobal_elems = get_sum(mesh->comm(), mesh->owned(dim));
  auto out = Write<Real>(mesh->nelems() * symm_ncomps(dim));
  auto f = LAMBDA(LO elem) {
    Few<LO, max_elems_per_patch> patch_elems;
    Int npatch_elems = 0;
    for (auto ev = elem * nverts_per_elem; ev < ((elem + 1) * nverts_per_elem);
         ++ev) {
      auto vert = elem_verts2vert[ev];
      for (auto ve = verts2elems.a2ab[vert]; ve < verts2elems.a2ab[vert + 1];
           ++ve) {
        auto patch_elem = verts2elems.ab2b[ve];
        add_unique(patch_elems, npatch_elems, patch_elem);
      }
    }
    Real patch_volume = 0;
    auto grad_sum = zero_vector<dim>();
    for (Int i = 0; i < npatch_elems; ++i) {
      auto patch_elem = patch_elems[i];
      auto gradient = get_vector<dim>(elem_gradients, patch_elem);
      auto volume = elems2volume[patch_elem];
      patch_volume += volume;
      grad_sum = grad_sum + (gradient * volume);
    }
    auto grad_avg = grad_sum / patch_volume;
    auto op_sum = zero_matrix<dim, dim>();
    for (Int i = 0; i < npatch_elems; ++i) {
      auto patch_elem = patch_elems[i];
      auto gradient = get_vector<dim>(elem_gradients, patch_elem);
      auto volume = elems2volume[patch_elem];
      auto grad_err = grad_avg - gradient;
      auto op = outer_product(grad_err, grad_err);
      op_sum = op_sum + (op * volume);
    }
    auto op_avg = op_sum / patch_volume;
    auto iso_volume =
        (dim == 3) ? (1.0 / (6.0 * sqrt(2.0))) : (sqrt(3.0) / 4.0);
    auto volume_factor = elems2volume[elem] / iso_volume;
    auto pullback_volume = patch_volume / volume_factor;
    auto a =
        square(error_bound) / (Real(dim) * nglobal_elems * pullback_volume);
    auto op_decomp = decompose_eigen(op_avg);
    auto g = op_decomp.l;
    auto gv = op_decomp.q;
    auto b = root<dim>(a);
    auto g_min = square(b / max_size);
    for (Int i = 0; i < dim; ++i) g[i] = max2(g[i], g_min);
    Matrix<dim, dim> r;
    for (Int i = 0; i < dim; ++i) r[i] = gv[dim - i - 1];
    Vector<dim> h;
    for (Int i = 0; i < dim; ++i) h[i] = b / sqrt(g[dim - i - 1]);
    auto m = compose_metric(r, h);
    set_symm(out, elem, m);
  };
  parallel_for(mesh->nelems(), f);
  return out;
}

Reals get_aniso_zz_metric(
    Mesh* mesh, Reals elem_gradients, Real error_bound, Real max_size) {
  if (mesh->dim() == 3) {
    return get_aniso_zz_metric_dim<3>(
        mesh, elem_gradients, error_bound, max_size);
  } else if (mesh->dim() == 2) {
    return get_aniso_zz_metric_dim<2>(
        mesh, elem_gradients, error_bound, max_size);
  } else {
    NORETURN(Reals());
  }
}

}  // end namespace Omega_h
