#include "Omega_h_file.hpp"
#include "Omega_h_build.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_vector.hpp"
#include "Omega_h_mesh.hpp"

#include "Omega_h_for.hpp"
#include "Omega_h_adj.hpp"

#include "SimPartitionedMesh.h"
#include "SimModel.h"
#include "SimUtil.h"

#include "SimDiscrete.h"

namespace Omega_h {

namespace meshsim {

int classId(pEntity e) {
  pGEntity g = EN_whatIn(e);
  assert(g);
  return GEN_tag(g);
}

void call_print(LOs a) {
  printf("\n");
  auto a_w = Write<LO> (a.size());
  auto r2w = OMEGA_H_LAMBDA(LO i) {
    a_w[i] = a[i];
  };
  parallel_for(a.size(), r2w);
  auto a_host = HostWrite<LO>(a_w);
  for (int i=0; i<a_host.size(); ++i) {
    printf(" %d,", a_host[i]);
  };
  printf("\n");
  printf("\n");
  return;
}

void read_internal(pMesh m, Mesh* mesh) {

  (void)mesh;

  const int numVtx = M_numVertices(m);
  const int numEdges = M_numEdges(m);
  const int numFaces = M_numFaces(m);
  const int numRegions = M_numRegions(m);

  std::vector<int> elem_vertices[4];
  std::vector<int> face_vertices[2];
  std::vector<int> edge_vertices[1];
  /* TODO:transfer classification info */

  Int max_dim;
  if (numRegions) {
    max_dim = 3;
  } else if (numFaces) {
    max_dim = 2;
  } else if (numEdges) {
    max_dim = 1;
  } else {
    Omega_h_fail("There were no Elements of dimension higher than zero!\n");
  }

  HostWrite<Real> host_coords(numVtx*max_dim);
  VIter vertices = M_vertexIter(m);
  pVertex vtx;
  LO v = 0;
  while ((vtx = (pVertex) VIter_next(vertices))) {
    double xyz[3];
    V_coord(vtx,xyz);
    if( max_dim < 3 && xyz[2] != 0 )
      Omega_h_fail("The z coordinate must be zero for a 2d mesh!\n");
    for(int j=0; j<max_dim; j++) {
      host_coords[v * max_dim + j] = xyz[j];
    }
    ++v;
  }
  VIter_delete(vertices);

  mesh->set_dim(max_dim);
  mesh->set_verts_type(numVtx);
  mesh->add_coords_mix(host_coords.write());

  edge_vertices[0].reserve(numEdges*2);
  EIter edges = M_edgeIter(m);
  pEdge edge;
  int count_edge = 0;
  while ((edge = (pEdge) EIter_next(edges))) {
    double xyz[3];
    count_edge += 1;
    for(int j=0; j<2; ++j) {
      vtx = E_vertex(edge,j);
      edge_vertices[0].push_back(EN_id(vtx));
      V_coord(vtx,xyz);
    }
  }
  EIter_delete(edges);
  
  HostWrite<LO> host_e2v(numEdges*2);
  for (Int i = 0; i < numEdges; ++i) {
    for (Int j = 0; j < 2; ++j) {
      host_e2v[i*2 + j] =
          edge_vertices[0][static_cast<std::size_t>(i*2 + j)];
    }
  }
  auto ev2v = Read<LO>(host_e2v.write());
  mesh->set_ents(Topo_type::edge, Topo_type::vertex, Adj(ev2v));

  FIter faces = M_faceIter(m);
  pFace face;
  int count_tri = 0;
  int count_quad = 0;
  while (face = (pFace) FIter_next(faces)) {
    if (F_numEdges(face) == 3) {
      count_tri += 1;
    }
    else if (F_numEdges(face) == 4) {
      count_quad += 1;
    }
    else {
      Omega_h_fail ("Face is neither tri nor quad \n");
    }
  }
  FIter_delete(faces);

  face_vertices[0].reserve(count_tri*3);
  face_vertices[1].reserve(count_quad*4);

  faces = M_faceIter(m);
  while (face = (pFace) FIter_next(faces)) {
    if (F_numEdges(face) == 3) {
      pVertex tri_vertex;
      pPList tri_vertices = F_vertices(face,1);
      assert (PList_size(tri_vertices) == 3);
      void *iter = 0;
      while (tri_vertex = (pVertex) PList_next(tri_vertices, &iter)) {
        face_vertices[0].push_back(EN_id(tri_vertex));
      }
      PList_delete(tri_vertices);
    }
    else if (F_numEdges(face) == 4) {
      pVertex quad_vertex;
      pPList quad_vertices = F_vertices(face,1);
      assert (PList_size(quad_vertices) == 4);
      void *iter = 0;
      while (quad_vertex = (pVertex) PList_next(quad_vertices, &iter)) {
        face_vertices[1].push_back(EN_id(quad_vertex));
      }
      PList_delete(quad_vertices);
    }
    else {
      Omega_h_fail ("Face is neither tri nor quad \n");
    }
  }
  FIter_delete(faces);
  auto edge2vert = mesh->get_adj(Topo_type::edge, Topo_type::vertex);
  auto vert2edge = mesh->ask_up(Topo_type::vertex, Topo_type::edge);
  HostWrite<LO> host_tri2verts(count_tri*3);
  for (Int i = 0; i < count_tri; ++i) {
    for (Int j = 0; j < 3; ++j) {
      host_tri2verts[i*3 + j] =
          face_vertices[0][static_cast<std::size_t>(i*3 + j)];
    }
  }
  auto tri2verts = Read<LO>(host_tri2verts.write());
  auto down = reflect_down(tri2verts, edge2vert.ab2b, vert2edge, Topo_type::triangle, Topo_type::edge);
  mesh->set_ents(Topo_type::triangle, Topo_type::edge, down);

  HostWrite<LO> host_quad2verts(count_quad*4);
  for (Int i = 0; i < count_quad; ++i) {
    for (Int j = 0; j < 4; ++j) {
      host_quad2verts[i*4 + j] =
          face_vertices[1][static_cast<std::size_t>(i*4 + j)];
    }
  }
  auto quad2verts = Read<LO>(host_quad2verts.write());
  down = reflect_down(quad2verts, edge2vert.ab2b, vert2edge, Topo_type::quadrilateral, Topo_type::edge);
  mesh->set_ents(Topo_type::quadrilateral, Topo_type::edge, down);

  RIter regions = M_regionIter(m);
  LO count_tet = 0;
  LO count_hex = 0;
  LO count_wedge = 0;
  LO count_pyramid = 0;
  pRegion rgn;
  while (rgn = (pRegion) RIter_next(regions)) {
    if (R_topoType(rgn) == Rtet) {
      count_tet += 1;
    }
    else if (R_topoType(rgn) == Rhex) {
      count_hex += 1;
    }
    else if (R_topoType(rgn) == Rwedge) {
      count_wedge += 1;
    }
    else if (R_topoType(rgn) == Rpyramid) {
      count_pyramid += 1;
    }
    else {
      Omega_h_fail("Region is not tet, hex, wedge, or pyramid \n");
    }
  }
  RIter_delete(regions);
  printf("tet=%d, hex=%d, wedge=%d, pyramid=%d\n",
         count_tet, count_hex, count_wedge, count_pyramid);

  elem_vertices[0].reserve(count_tet*4);
  elem_vertices[1].reserve(count_hex*8);
  elem_vertices[2].reserve(count_wedge*6);
  elem_vertices[3].reserve(count_pyramid*5);

  regions = M_regionIter(m);
  while ((rgn = (pRegion) RIter_next(regions))) {
    if (R_topoType(rgn) == Rtet) {
      pVertex vert;
      pPList verts = R_vertices(rgn,1);
      assert (PList_size(verts) == 4);
      void *iter = 0;
      while (vert = (pVertex) PList_next(verts, &iter)) {
        elem_vertices[0].push_back(EN_id(vert));
      }
      PList_delete(verts);
    }
    else if (R_topoType(rgn) == Rhex) {
      pVertex vert;
      pPList verts = R_vertices(rgn,1);
      assert (PList_size(verts) == 8);
      void *iter = 0;
      while (vert = (pVertex) PList_next(verts, &iter)) {
        elem_vertices[1].push_back(EN_id(vert));
      }
      PList_delete(verts);
    }
    else if (R_topoType(rgn) == Rwedge) {
      pVertex vert;
      pPList verts = R_vertices(rgn,1);
      assert (PList_size(verts) == 6);
      void *iter = 0;
      while (vert = (pVertex) PList_next(verts, &iter)) {
        elem_vertices[2].push_back(EN_id(vert));
      }
      PList_delete(verts);
    }
    else if (R_topoType(rgn) == Rpyramid) {
      pVertex vert;
      pPList verts = R_vertices(rgn,1);
      assert (PList_size(verts) == 5);
      void *iter = 0;
      while (vert = (pVertex) PList_next(verts, &iter)) {
        elem_vertices[3].push_back(EN_id(vert));
      }
      PList_delete(verts);
    }
    else {
      Omega_h_fail ("Region is not tet, hex, wedge, or pyramid \n");
    }
  }
  RIter_delete(regions);

  auto tri2vert = mesh->ask_down(Topo_type::triangle, Topo_type::vertex);
  auto vert2tri = mesh->ask_up(Topo_type::vertex, Topo_type::triangle);
  auto quad2vert = mesh->ask_down(Topo_type::quadrilateral, Topo_type::vertex);
  auto vert2quad = mesh->ask_up(Topo_type::vertex, Topo_type::quadrilateral);

  HostWrite<LO> host_tet2verts(count_tet*4);
  for (Int i = 0; i < count_tet; ++i) {
    for (Int j = 0; j < 4; ++j) {
      host_tet2verts[i*4 + j] =
          elem_vertices[0][static_cast<std::size_t>(i*4 + j)];
    }
  }
  auto tet2verts = Read<LO>(host_tet2verts.write());
  down = reflect_down(tet2verts, tri2vert.ab2b, vert2tri, Topo_type::tetrahedron, Topo_type::triangle);
  mesh->set_ents(Topo_type::tetrahedron, Topo_type::triangle, down);

  HostWrite<LO> host_hex2verts(count_hex*8);
  for (Int i = 0; i < count_hex; ++i) {
    for (Int j = 0; j < 8; ++j) {
      host_hex2verts[i*8 + j] =
          elem_vertices[1][static_cast<std::size_t>(i*8 + j)];
    }
  }
  auto hex2verts = Read<LO>(host_hex2verts.write());
  down = reflect_down(hex2verts, quad2vert.ab2b, vert2quad, Topo_type::hexahedron, Topo_type::quadrilateral);
  mesh->set_ents(Topo_type::hexahedron, Topo_type::quadrilateral, down);

  HostWrite<LO> host_wedge2verts(count_wedge*6);
  for (Int i = 0; i < count_wedge; ++i) {
    for (Int j = 0; j < 6; ++j) {
      host_wedge2verts[i*6 + j] =
          elem_vertices[2][static_cast<std::size_t>(i*6 + j)];
    }
  }
  auto wedge2verts = Read<LO>(host_wedge2verts.write());
  down = reflect_down(wedge2verts, quad2vert.ab2b, vert2quad, Topo_type::wedge, Topo_type::quadrilateral);
  mesh->set_ents(Topo_type::wedge, Topo_type::quadrilateral, down);

  down = reflect_down(wedge2verts, tri2vert.ab2b, vert2tri, Topo_type::wedge, Topo_type::triangle);
  mesh->set_ents(Topo_type::wedge, Topo_type::triangle, down);

  HostWrite<LO> host_pyramid2verts(count_pyramid*5);
  for (Int i = 0; i < count_pyramid; ++i) {
    for (Int j = 0; j < 5; ++j) {
      host_pyramid2verts[i*5 + j] =
          elem_vertices[3][static_cast<std::size_t>(i*5 + j)];
    }
  }
  auto pyramid2verts = Read<LO>(host_pyramid2verts.write());
  down = reflect_down(pyramid2verts, tri2vert.ab2b, vert2tri, Topo_type::pyramid, Topo_type::triangle);
  mesh->set_ents(Topo_type::pyramid, Topo_type::triangle, down);

  down = reflect_down(pyramid2verts, quad2vert.ab2b, vert2quad, Topo_type::pyramid, Topo_type::quadrilateral);
  mesh->set_ents(Topo_type::pyramid, Topo_type::quadrilateral, down);

}

Mesh read(filesystem::path const& mesh_fname, filesystem::path const& mdl_fname,
    CommPtr comm) {
  SimPartitionedMesh_start(NULL,NULL);
  SimModel_start();
  Sim_readLicenseFile(NULL);
  SimDiscrete_start(0);
  pNativeModel nm = NULL;
  pProgress p = NULL;
  pGModel g = GM_load(mdl_fname.c_str(), nm, p);
  pMesh m = M_load(mesh_fname.c_str(), g, p);
  auto mesh = Mesh(comm->library());
  meshsim::read_internal(m, &mesh);
  M_release(m);
  GM_release(g);
  SimDiscrete_stop(0);
  SimModel_stop();
  return mesh;
}

}  // namespace meshsim

}  // end namespace Omega_h
