setZ0
Aplication for setting nonUniform z0 (C++). 
It is based on foamDictionary, an application for managing dictionary entries.
Author: Christina Fratzeskou

Tested on OpenFoam 7.

Requirements:
- An installation of OpenFoam 7
- An obj file for the input geometry, with landcover classification
- A .mtl file with landcover classification

Installation:
- Download z0tool
- Compile using wmake
- Copy setZ0Dict to the /constant folder of your case,
  and modify it (Make sure landcover classification match 
  the one in the .obj and .mtl files).

Example usage:
Setting z0 to a ground patch:
setZ0 0/nut -entry boundaryField.terrain.z0 -setZ0Ground terrain 
-writeCoords -writeZ0

Set z0 to inlet patch:
setZ0 0/epsilon -entry boundaryField.inlet.z0 -setZ0Inlet inlet 
-writeCoords -writeZ0
