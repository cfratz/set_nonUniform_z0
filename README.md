setZ0:

Aplication for setting nonUniform z0 to a patch(C++). 
Author: Christina Fratzeskou

Tested on OpenFoam v7 on Ubuntu 18.04LTS on Windows10.

Requirements:
- An installation of OpenFoam v7

Installation:
- Download setZ0.C, Make folder
- Add the Make and setZ0 file to a directory
- Compile using wmake

Case preparation:
- Copy setZ0Dict to the /constant folder of your case,
  and modify it (Make sure landcover classification match 
  the one in the .obj and .mtl files).
 - Prepare boundary and initial conditions
 - Mesh generation (Suggested OpenFOAM mesh tools: blockMesh, surfaceFeatures, snappyHexMesh)
 
Input:
- An obj file for the input geometry, with landcover classification
- A .mtl file with landcover classification

Example usage:
  Setting z0 to a ground patch:
  setZ0 0/nut -entry boundaryField.terrain.z0 -setZ0Ground terrain [optional] -setZ0NoGeom 0.5 -writeCoords -writeZ0 -exportToVtk

  Set z0 to inlet patch:
  setZ0 0/epsilon -entry boundaryField.inlet.z0 -setZ0Inlet inlet [optional] -setParams -writeCoords -writeZ0 -exportToVtk
