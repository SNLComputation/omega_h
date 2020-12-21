#include <iostream>
#include <string>

#include "Omega_h_file.hpp"
#include "Omega_h_mesh.hpp"

#include "Omega_h_element.hpp"
#include "Omega_h_array_ops.hpp"
using namespace Omega_h;

void test_degree() {
  OMEGA_H_CHECK(element_degree(Topo_type::edge, Topo_type::vertex) == 2);
  OMEGA_H_CHECK(element_degree(Topo_type::tetrahedron, Topo_type::vertex) == 4);
  OMEGA_H_CHECK(element_degree(Topo_type::hexahedron, Topo_type::quadrilateral) == 6);
  OMEGA_H_CHECK(element_degree(Topo_type::wedge, Topo_type::triangle) == 2);
  OMEGA_H_CHECK(element_degree(Topo_type::pyramid, Topo_type::edge) == 8);
  OMEGA_H_CHECK(element_degree(Topo_type::pyramid, Topo_type::quadrilateral) == 1);
}

void test_tags(Mesh* mesh) {
  auto num_wedge = mesh->nwedges();
  mesh->add_tag<LO>(Topo_type::wedge, "gravity", 1);
  mesh->set_tag<LO>(Topo_type::wedge, "gravity", LOs(num_wedge,10));
  auto test_tag1 = mesh->get_array<LO>(Topo_type::wedge, "gravity");
  OMEGA_H_CHECK(test_tag1 == LOs(num_wedge, 10));
  mesh->remove_tag(Topo_type::wedge, "gravity");

  auto num_pyram = mesh->npyrams();
  mesh->add_tag<Real>(Topo_type::pyramid, "density", 1, Reals(num_pyram, 0.0005));
  auto test_tag2 = mesh->get_array<Real>(Topo_type::pyramid, "density");
  OMEGA_H_CHECK(test_tag2 == Reals(num_pyram, 0.0005));
  mesh->remove_tag(Topo_type::pyramid, "density");
}

void test_adjs(Mesh* mesh) {

  //get number of entities
  auto num_vertex = mesh->nverts_mix();
  auto num_edge = mesh->nedges_mix();
  auto num_tri = mesh->ntris();
  auto num_quad = mesh->nquads();
  auto num_tet = mesh->ntets();
  auto num_hex = mesh->nhexs();
  auto num_wedge = mesh->nwedges();
  auto num_pyramid = mesh->npyrams();

  //1 lvl dwn queries
  auto edge2vert = mesh->get_adj(Topo_type::edge, Topo_type::vertex);
  auto tri2edge = mesh->get_adj(Topo_type::triangle, Topo_type::edge);
  auto quad2edge = mesh->get_adj(Topo_type::quadrilateral, Topo_type::edge);
  auto tet2tri = mesh->get_adj(Topo_type::tetrahedron, Topo_type::triangle);
  auto hex2quad = mesh->get_adj(Topo_type::hexahedron, Topo_type::quadrilateral);
  auto wedge2tria = mesh->get_adj(Topo_type::wedge, Topo_type::triangle);
  auto wedge2quadr = mesh->get_adj(Topo_type::wedge, Topo_type::quadrilateral);
  auto pyramid2tria = mesh->get_adj(Topo_type::pyramid, Topo_type::triangle);
  auto pyramid2quadr = mesh->get_adj(Topo_type::pyramid, Topo_type::quadrilateral);

  //1 lvl inversion queries
  auto vert2edge = mesh->ask_up(Topo_type::vertex, Topo_type::edge);
  auto edge2tri = mesh->ask_up(Topo_type::edge, Topo_type::triangle);
  auto edge2quad = mesh->ask_up(Topo_type::edge, Topo_type::quadrilateral);
  auto tri2tet = mesh->ask_up(Topo_type::triangle, Topo_type::tetrahedron);
  auto quad2hex = mesh->ask_up(Topo_type::quadrilateral, Topo_type::hexahedron);
  auto tri2wedge = mesh->ask_up(Topo_type::triangle, Topo_type::wedge);
  auto quad2wedge = mesh->ask_up(Topo_type::quadrilateral, Topo_type::wedge);
  auto tri2pyramid = mesh->ask_up(Topo_type::triangle, Topo_type::pyramid);
  auto quad2pyramid = mesh->ask_up(Topo_type::quadrilateral, Topo_type::pyramid);

  //1-lvl inversion assertions
  OMEGA_H_CHECK(vert2edge.ab2b.size() == 2*num_edge);
  OMEGA_H_CHECK(vert2edge.a2ab.size() == num_vertex+1);
  OMEGA_H_CHECK(edge2tri.ab2b.size() == 3*num_tri);
  OMEGA_H_CHECK(edge2tri.a2ab.size() == num_edge+1);
  OMEGA_H_CHECK(edge2quad.ab2b.size() == 4*num_quad);
  OMEGA_H_CHECK(edge2quad.a2ab.size() == num_edge+1);
  OMEGA_H_CHECK(tri2tet.ab2b.size() == 4*num_tet);
  OMEGA_H_CHECK(tri2tet.a2ab.size() == num_tri+1);
  OMEGA_H_CHECK(quad2hex.ab2b.size() == num_hex*6);
  OMEGA_H_CHECK(quad2hex.a2ab.size() == num_quad+1);
  OMEGA_H_CHECK(tri2wedge.ab2b.size() == 2*num_wedge);
  OMEGA_H_CHECK(tri2wedge.a2ab.size() == num_tri+1);
  OMEGA_H_CHECK(quad2wedge.ab2b.size() == num_wedge*3);
  OMEGA_H_CHECK(quad2wedge.a2ab.size() == num_quad+1);
  OMEGA_H_CHECK(tri2pyramid.ab2b.size() == 4*num_pyramid);
  OMEGA_H_CHECK(tri2pyramid.a2ab.size() == num_tri+1);
  OMEGA_H_CHECK(quad2pyramid.ab2b.size() == num_pyramid);
  OMEGA_H_CHECK(quad2pyramid.a2ab.size() == num_quad+1);

  //transit tests
  auto tri2vert = mesh->ask_down(Topo_type::triangle, Topo_type::vertex);
  auto quad2vert = mesh->ask_down(Topo_type::quadrilateral, Topo_type::vertex);
  auto tet2edge = mesh->ask_down(Topo_type::tetrahedron, Topo_type::edge);
  auto tet2vtx = mesh->ask_down(Topo_type::tetrahedron, Topo_type::vertex);
  auto hex2edge = mesh->ask_down(Topo_type::hexahedron, Topo_type::edge);
  auto hex2vtx = mesh->ask_down(Topo_type::hexahedron, Topo_type::vertex);
  auto wedge2edge = mesh->ask_down(Topo_type::wedge, Topo_type::edge);
  auto wedge2vtx = mesh->ask_down(Topo_type::wedge, Topo_type::vertex);
  auto pyram2edge = mesh->ask_down(Topo_type::pyramid, Topo_type::edge);
  auto pyram2vtx = mesh->ask_down(Topo_type::pyramid, Topo_type::vertex);

  //upward
  auto vert2tri = mesh->ask_up(Topo_type::vertex, Topo_type::triangle);
  auto vert2quad = mesh->ask_up(Topo_type::vertex, Topo_type::quadrilateral);
  auto vert2tet = mesh->ask_up(Topo_type::vertex, Topo_type::tetrahedron);
  auto vert2hex = mesh->ask_up(Topo_type::vertex, Topo_type::hexahedron);
  auto vert2wedge = mesh->ask_up(Topo_type::vertex, Topo_type::wedge);
  auto vert2pyramid = mesh->ask_up(Topo_type::vertex, Topo_type::pyramid);
  auto edge2tet = mesh->ask_up(Topo_type::edge, Topo_type::tetrahedron);
  auto edge2hex = mesh->ask_up(Topo_type::edge, Topo_type::hexahedron);
  auto edge2wedge = mesh->ask_up(Topo_type::edge, Topo_type::wedge);
  auto edge2pyramid = mesh->ask_up(Topo_type::edge, Topo_type::pyramid);

  //transit assertions
  OMEGA_H_CHECK(tri2vert.ab2b.size() == num_tri*3);
  OMEGA_H_CHECK(quad2vert.ab2b.size() == num_quad*4);
  OMEGA_H_CHECK(tet2edge.ab2b.size() == num_tet*6);
  OMEGA_H_CHECK(tet2vtx.ab2b.size() == num_tet*4);
  OMEGA_H_CHECK(hex2edge.ab2b.size() == num_hex*12);
  OMEGA_H_CHECK(hex2vtx.ab2b.size() == num_hex*8);
  OMEGA_H_CHECK(wedge2edge.ab2b.size() == num_wedge*9);
  OMEGA_H_CHECK(wedge2vtx.ab2b.size() == num_wedge*6);
  OMEGA_H_CHECK(pyram2edge.ab2b.size() == num_pyramid*8);
  OMEGA_H_CHECK(pyram2vtx.ab2b.size() == num_pyramid*5);

  //check values for example meshes
  if (num_vertex == 12) {
    //4 elem mesh
    OMEGA_H_CHECK(tet2vtx.ab2b == LOs({0, 1, 8, 2}));
    OMEGA_H_CHECK(hex2vtx.ab2b == LOs({4, 5, 9, 11, 7, 6, 1, 0}));
    OMEGA_H_CHECK(wedge2vtx.ab2b == LOs({11, 9, 10, 0, 1, 8}));
    OMEGA_H_CHECK(pyram2vtx.ab2b == LOs({9, 10, 8, 1, 3}));
    OMEGA_H_CHECK(tet2edge.ab2b == LOs({0, 22, 20, 1, 3, 6}));
    OMEGA_H_CHECK(hex2edge.ab2b == LOs({13, 9, 17, 12, 14, 15, 23, 21, 16, 5, 0, 2}));
    OMEGA_H_CHECK(wedge2edge.ab2b == LOs({17, 18, 11, 21, 23, 19, 0, 22, 20}));
    OMEGA_H_CHECK(pyram2edge.ab2b == LOs({18, 19, 22, 23, 8, 10, 7, 4}));
  }
  else if (num_vertex == 9) {
    //pyram on hex
    OMEGA_H_CHECK(hex2edge.ab2b == LOs({0, 4, 7, 1, 2, 5, 8, 10, 12, 13, 14, 15}));
    OMEGA_H_CHECK(hex2vtx.ab2b == LOs({0, 1, 2, 3, 5, 6, 7, 8}));
    OMEGA_H_CHECK(pyram2edge.ab2b == LOs({12, 13, 14, 15, 3, 6, 9, 11}));
    OMEGA_H_CHECK(pyram2vtx.ab2b == LOs({5, 6, 7, 8, 4}));
  }
  else if (num_vertex == 7) {
    //tet on wedge
    OMEGA_H_CHECK(tet2edge.ab2b == LOs({9, 10, 11, 3, 6, 8})); 
    OMEGA_H_CHECK(tet2vtx.ab2b == LOs({4, 5, 6, 3})); 
    OMEGA_H_CHECK(wedge2edge.ab2b == LOs({0, 4, 1, 2, 5, 7, 9, 10, 11}));
    OMEGA_H_CHECK(wedge2vtx.ab2b == LOs({0, 1, 2, 4, 5, 6})); 
  }
}

void test_finer_meshes(CommPtr comm) {

  //tet=4609, hex=249, wedge=262, pyramid=313
  std::string mesh_in = "/users/joshia5/simmodeler/localconcave_turorial_mixedvol_geomsim3-case1.sms";
  std::string model_in = "/users/joshia5/simmodeler/localconcave_turorial_mixedvol_geomsim3.smd";
  auto mesh = meshsim::read(mesh_in, model_in, comm);
  test_adjs(&mesh);
  test_tags(&mesh);
  vtk::write_vtu("/users/joshia5/simmodeler/geomsim3_4type.vtu", &mesh, Topo_type::pyramid);

  //tet=1246, hex=8, wedge=0, pyramid=229
  mesh_in = "/users/joshia5/simmodeler/localconcave_turorial_mixedvol_geomsim2-case1.sms";
  model_in = "/users/joshia5/simmodeler/localconcave_turorial_mixedvol_geomsim2.smd";
  mesh = meshsim::read(mesh_in, model_in, comm);
  test_adjs(&mesh);
  test_tags(&mesh);

  //tet=3274, hex=68, wedge=0, pyramid=171
  mesh_in = "/users/joshia5/simmodeler/localconcave_turorial_mixedvol_geomsim-case1.sms";
  model_in = "/users/joshia5/simmodeler/localconcave_turorial_mixedvol_geomsim.smd";
  mesh = meshsim::read(mesh_in, model_in, comm);
  test_adjs(&mesh);
  test_tags(&mesh);
  vtk::write_vtu("/users/joshia5/simmodeler/geomsim_tet_he_py.vtu", &mesh, Topo_type::pyramid);

  //tet=4437, hex=0, wedge=0, pyramid=0
  mesh_in = "/users/joshia5/simmodeler/localconcave_turorial-case1_v7.sms";
  model_in = "/users/joshia5/simmodeler/localconcave_turorial_geomsim.smd";
  mesh = meshsim::read(mesh_in, model_in, comm);
  test_adjs(&mesh);
  test_tags(&mesh);
  vtk::write_vtu("/users/joshia5/simmodeler/geomsim_tet.vtu", &mesh, Topo_type::pyramid);

}

int main(int argc, char** argv) {
  test_degree();

  auto lib = Library(&argc, &argv);
  auto comm = lib.world();

  /* Generate these meshes using the mixed_writeMesh ctest */
  std::string mesh_in = "/users/joshia5/simmodeler/Example_hex.sms";
  std::string model_in = "/users/joshia5/simmodeler/Example_hex.smd";
  auto mesh = meshsim::read(mesh_in, model_in, comm);
  test_adjs(&mesh);
  test_tags(&mesh);
  vtk::write_vtu("/users/joshia5/simmodeler/hex.vtu", &mesh, Topo_type::pyramid);

  mesh_in = "/users/joshia5/simmodeler/Example_wedge.sms";
  model_in = "/users/joshia5/simmodeler/Example_wedge.smd";
  mesh = meshsim::read(mesh_in, model_in, comm);
  test_adjs(&mesh);
  test_tags(&mesh);
  vtk::write_vtu("/users/joshia5/simmodeler/wedge.vtu", &mesh, Topo_type::pyramid);

  mesh_in = "/users/joshia5/simmodeler/Example_pym.sms";
  model_in = "/users/joshia5/simmodeler/Example_pym.smd";
  mesh = meshsim::read(mesh_in, model_in, comm);
  test_adjs(&mesh);
  test_tags(&mesh);
  vtk::write_vtu("/users/joshia5/simmodeler/pym.vtu", &mesh, Topo_type::pyramid);

  mesh_in = "/users/joshia5/simmodeler/Example_tet_wedge.sms";
  model_in = "/users/joshia5/simmodeler/Example_tet_wedge.smd";
  mesh = meshsim::read(mesh_in, model_in, comm);
  test_adjs(&mesh);
  test_tags(&mesh);
  vtk::write_vtu("/users/joshia5/simmodeler/tet_wedge.vtu", &mesh, Topo_type::pyramid);

  mesh_in = "/users/joshia5/simmodeler/Example_pym_hex.sms";
  model_in = "/users/joshia5/simmodeler/Example_pym_hex.smd";
  mesh = meshsim::read(mesh_in, model_in, comm);
  test_adjs(&mesh);
  test_tags(&mesh);
  vtk::write_vtu("/users/joshia5/simmodeler/pym_hex.vtu", &mesh, Topo_type::pyramid);

  mesh_in = "/users/joshia5/simmodeler/Example_allType.sms";
  model_in = "/users/joshia5/simmodeler/Example_allType.smd";
  mesh = meshsim::read(mesh_in, model_in, comm);
  test_adjs(&mesh);
  test_tags(&mesh);
  vtk::write_vtu("/users/joshia5/simmodeler/4elems.vtu", &mesh, Topo_type::pyramid);

  /* call if the finer test meshes are available */
  /* test_finer_meshes(comm); */

  return 0;
}
