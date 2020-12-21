#include "MeshSim.h"
#include "SimDiscrete.h"
#include "SimMessages.h"
#include "SimError.h"
#include "SimErrorCodes.h"
#include "SimMeshingErrorCodes.h"
#include "SimDiscreteErrorCodes.h"
#include <iostream>
using namespace std;

#include "Omega_h_comm.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_file.hpp"
using namespace Omega_h;

void messageHandler(int type, const char *msg)
{
  switch (type) {
  case Sim_InfoMsg:
    cout<<"Info: "<<msg<<endl;
    break;
  case Sim_DebugMsg:
    cout<<"Debug: "<<msg<<endl;
    break;
  case Sim_WarningMsg:
    cout<<"Warning: "<<msg<<endl;
    break;
  case Sim_ErrorMsg:
    cout<<"Error: "<<msg<<endl;
    break;
  }
  return;
}

void finalize_write(int numVerts, const double *coords, int numElems,
                    int *elementType, int *elementData, pVertex *vReturn,
                    pEntity *eReturn, const char *mesh_path,
                    const char *model_path) {
  Sim_logOn("importData1.log");
  MS_init();
  Sim_readLicenseFile(0);
  pMesh meshtest;
  SimDiscrete_start(0);
  Sim_setMessageHandler(messageHandler);
  pProgress progress = Progress_new();
  Progress_setDefaultCallback(progress);

  meshtest = M_new(0,0);
  if(M_importFromData(meshtest, numVerts, coords, numElems,
                      elementType, elementData, vReturn, eReturn,
                      progress)) {
    cerr<<"Error importing mesh data"<<endl;
    M_release(meshtest);
  }

  pDiscreteModel modeltest = DM_createFromMesh(meshtest, 0, progress);
  if(!modeltest) {
    cerr<<"Error creating Discrete model from mesh"<<endl;
    M_release(meshtest);
  }

  DM_findEdgesByFaceNormals(modeltest, 0, progress);
  DM_eliminateDanglingEdges(modeltest, progress);
  if(DM_completeTopology(modeltest, progress)) {
    cerr<<"Error completing Discrete model topology"<<endl;
    M_release(meshtest);
    GM_release(modeltest);
  }

  GM_write(modeltest, model_path, 0, progress);
  M_write(meshtest, mesh_path, 0, progress);

  M_release(meshtest);
  GM_release(modeltest);
  Progress_delete(progress);
  SimDiscrete_stop(0);
  Sim_unregisterAllKeys();
  MS_exit();
  Sim_logOff();
}

void hex_gen(CommPtr comm) {
  auto mesh = Mesh(comm->library());
  try {
    int numVerts = 8;
    double coords[8*3] =  {
                           0.000, 0.000, 0.000,
                           1.000, 0.000, 0.000,
                           1.000, 1.000, 0.000,
                           0.000, 1.000, 0.000,
                           0.000, 0.000, 1.000,
                           1.000, 0.000, 1.000,
                           1.000, 1.000, 1.000,
                           0.000, 1.000, 1.000
                          };
    int numElems = 1;
    int elementType[1] = {13};
    int elementData[8+5] = {0,1,2,3,4,5,6,7};
    pVertex vReturn[8];
    pEntity eReturn[1];
    const char *mesh_path = "/users/joshia5/simmodeler/Example_hex.sms";
    const char *model_path = "/users/joshia5/simmodeler/Example_hex.smd";
    finalize_write(numVerts, coords, numElems, elementType, elementData,
                   vReturn, eReturn, mesh_path, model_path);
  } catch (pSimError err) {
    cerr<<"SimModSuite error caught:"<<endl;
    cerr<<"  Error code: "<<SimError_code(err)<<endl;
    cerr<<"  Error string: "<<SimError_toString(err)<<endl;
    SimError_delete(err);
  } catch (...) {
    cerr<<"Unhandled exception caught"<<endl;
  }
}

void wedge_gen(CommPtr comm) {
  auto mesh = Mesh(comm->library());
  try {
    int numVerts = 6;
    double coords[6*3] =  {0.000, 0.000, 0.000,
                           1.000, 0.000, 0.000,
                           0.500, 0.866, 0.000,
                           0.000, 0.000, 1.000,
                           1.000, 0.000, 1.000,
                           0.500, 0.866, 1.000,
                          };
    int numElems = 1;
    int elementType[1] = {12};
    int elementData[6] = {0,1,2,3,4,5};
    pVertex vReturn[6];
    pEntity eReturn[1];
    const char *mesh_path = "/users/joshia5/simmodeler/Example_wedge.sms";
    const char *model_path = "/users/joshia5/simmodeler/Example_wedge.smd";
    finalize_write(numVerts, coords, numElems, elementType, elementData,
                   vReturn, eReturn, mesh_path, model_path);
  } catch (pSimError err) {
    cerr<<"SimModSuite error caught:"<<endl;
    cerr<<"  Error code: "<<SimError_code(err)<<endl;
    cerr<<"  Error string: "<<SimError_toString(err)<<endl;
    SimError_delete(err);
  } catch (...) {
    cerr<<"Unhandled exception caught"<<endl;
  }
}

void pyramid_gen(CommPtr comm) {
  auto mesh = Mesh(comm->library());
  try {
    int numVerts = 5;
    double coords[5*3] =  {
                           0.000, 0.000, 1.000,
                           1.000, 0.000, 1.000,
                           1.000, 1.000, 1.000,
                           0.000, 1.000, 1.000,
                           0.500, 0.500, 2.000
                          };
    int numElems = 1;
    int elementType[1] = {11};
    int elementData[5] = {0,1,2,3,4};
    pVertex vReturn[5];
    pEntity eReturn[1];
    const char *mesh_path = "/users/joshia5/simmodeler/Example_pym.sms";
    const char *model_path = "/users/joshia5/simmodeler/Example_pym.smd";
    finalize_write(numVerts, coords, numElems, elementType, elementData,
                   vReturn, eReturn, mesh_path, model_path);
  } catch (pSimError err) {
    cerr<<"SimModSuite error caught:"<<endl;
    cerr<<"  Error code: "<<SimError_code(err)<<endl;
    cerr<<"  Error string: "<<SimError_toString(err)<<endl;
    SimError_delete(err);
  } catch (...) {
    cerr<<"Unhandled exception caught"<<endl;
  }
}

void tetOnWedge_gen(CommPtr comm) {
  auto mesh = Mesh(comm->library());
  try {
    int numVerts = 7;
    double coords[7*3] =  {0.000, 0.000, 0.000,
                           1.000, 0.000, 0.000,
                           0.500, 0.866, 0.000,
                           0.000, 0.000, 1.000,
                           1.000, 0.000, 1.000,
                           0.500, 0.866, 1.000,
                           0.500, 0.289, 1.866,
                          };
    int numElems = 2;
    int elementType[2] = {12, 10};
    int elementData[6+4] = {0,1,2,3,4,5,3,4,5,6};
    pVertex vReturn[7];
    pEntity eReturn[2];
    const char *mesh_path = "/users/joshia5/simmodeler/Example_tet_wedge.sms";
    const char *model_path = "/users/joshia5/simmodeler/Example_tet_wedge.smd";
    finalize_write(numVerts, coords, numElems, elementType, elementData,
                   vReturn, eReturn, mesh_path, model_path);
  } catch (pSimError err) {
    cerr<<"SimModSuite error caught:"<<endl;
    cerr<<"  Error code: "<<SimError_code(err)<<endl;
    cerr<<"  Error string: "<<SimError_toString(err)<<endl;
    SimError_delete(err);
  } catch (...) {
    cerr<<"Unhandled exception caught"<<endl;
  }
}

void pymOnHex_gen(CommPtr comm) {
  auto mesh = Mesh(comm->library());
  try {
    int numVerts = 9;
    double coords[9*3] =  {
                           0.000, 0.000, 0.000,
                           1.000, 0.000, 0.000,
                           1.000, 1.000, 0.000,
                           0.000, 1.000, 0.000,
                           0.000, 0.000, 1.000,
                           1.000, 0.000, 1.000,
                           1.000, 1.000, 1.000,
                           0.000, 1.000, 1.000,
                           0.500, 0.500, 2.000
                          };
    int numElems = 2;
    int elementType[2] = {13, 11};
    int elementData[8+5] = {0,1,2,3,4,5,6,7, 4,5,6,7,8};
    pVertex vReturn[9];
    pEntity eReturn[2];
    const char *mesh_path = "/users/joshia5/simmodeler/Example_pym_hex.sms";
    const char *model_path = "/users/joshia5/simmodeler/Example_pym_hex.smd";
    finalize_write(numVerts, coords, numElems, elementType, elementData,
                   vReturn, eReturn, mesh_path, model_path);
  } catch (pSimError err) {
    cerr<<"SimModSuite error caught:"<<endl;
    cerr<<"  Error code: "<<SimError_code(err)<<endl;
    cerr<<"  Error string: "<<SimError_toString(err)<<endl;
    SimError_delete(err);
  } catch (...) {
    cerr<<"Unhandled exception caught"<<endl;
  }
}

void allType_gen(CommPtr comm) {
  auto mesh = Mesh(comm->library());
  try {
    int numVerts = 12;
    double coords[12*3] = {0.000, 0.000, 0.000,
                           1.000, 0.000, 0.000,
                           0.500, 0.866, 0.000,
                           0.000, 0.000, 1.000, 
                           1.000, 0.000, 1.000,
                           0.500, 0.866, 1.000,
                           0.500, 0.289, 1.866,
                           1.000, 1.500, 0.500,
                           0.000, -1.000, 0.000,
                           1.000, -1.000, 0.000,
                           0.000, -1.000, 1.000,
                           1.000, -1.000, 1.000
                          };
    int numElems = 4;
    int elementType[4] = {12, 10, 11, 13};
    int elementData[6+4+5+8] = {0,1,2,3,4,5, 3,4,5,6, 5,4,1,2,7,
8,9,1,0,10,11,4,3};
    pVertex vReturn[12]; 
    pEntity eReturn[4]; 
    const char *mesh_path = "/users/joshia5/simmodeler/Example_allType.sms";
    const char *model_path = "/users/joshia5/simmodeler/Example_allType.smd";
    finalize_write(numVerts, coords, numElems, elementType, elementData,
                   vReturn, eReturn, mesh_path, model_path);
  } catch (pSimError err) {
    cerr<<"SimModSuite error caught:"<<endl;
    cerr<<"  Error code: "<<SimError_code(err)<<endl;
    cerr<<"  Error string: "<<SimError_toString(err)<<endl;
    SimError_delete(err);
  } catch (...) {
    cerr<<"Unhandled exception caught"<<endl;
  }
}

int main(int argc, char *argv[]) {
  auto lib = Library(&argc, &argv);
  auto comm = lib.world();

  hex_gen(comm);
  wedge_gen(comm);
  pyramid_gen(comm);
  tetOnWedge_gen(comm);
  pymOnHex_gen(comm);
  allType_gen(comm);

  return 0;
}
