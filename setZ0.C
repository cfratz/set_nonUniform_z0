/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    setNonUniformZ0  based on foamDictionary found in link: 
    https://cpp.openfoam.org/v7/foamDictionary_8C_source.html 

Description
    Assigns roughness length(i.e. z0) values to patches with nonuniform z0 values.
    It can also be used to manipulate dictionaries.

Usage
    \b setZ0 [OPTION] dictionary
      
      From foamDictionary:

      - \par -entry \<name\>
        Selects an entry

      - \par -keywords \<name\>
        Prints the keywords (of the selected entry or of the top level if
        no entry was selected

      - \par -add \<value\>
        Adds the entry (should not exist yet)

      - \par -set \<value\>
        Adds or replaces the entry

      - \par -remove
        Remove the selected entry

      - \par -expand
        Read the specified dictionary file, expand the macros etc. and write
        the resulting dictionary to standard output.

      - \par -includes
        List the \c #include and \c #includeIfPresent files to standard output

      - \par -disableFunctionEntries
        Do not expand macros or directives (#include etc)
      
      - \par -writeCoords
        Writes face centers of input ground patch to .xyz file

      - \par -writeZ0
        Writes z0 list of input patch to .txt file

      - \par -setZ0Ground \<patchName\>
        Sets or replaces z0 to the specified ground patch.
	This is meant for cases where the roughness length(z0) values for the specified
        patch are non uniform. For this the input triangulated geometry is required to 
	have added semantics (i.e obj + uml).

      - \par -setZ0Inlet \<patchName\>
        Sets or replaces z0 value to the specified inlet patch.
        This is meant for cases where the roughness length(z0) values for the specified
	patch are non uniform.For this the input triangulated geometry is required to 
	have added semantics (i.e obj + uml).

      -\par -
    Example usage:
      - Change simulation to run for one timestep only:
        \verbatim
          foamDictionary system/controlDict -entry stopAt -set writeNow
        \endverbatim

      - Change solver:
        \verbatim
           foamDictionary system/fvSolution -entry solvers.p.solver -set PCG
        \endverbatim

      - Print bc type:
        \verbatim
           foamDictionary 0/U -entry boundaryField.movingWall.type
        \endverbatim

      - Change bc parameter:
        \verbatim
           foamDictionary 0/U -entry boundaryField.movingWall.value \
             -set "uniform (2 0 0)"
        \endverbatim

      - Change whole bc type:
        \verbatim
          foamDictionary 0/U -entry boundaryField.movingWall \
            -set "{type uniformFixedValue; uniformValue (2 0 0);}"
        \endverbatim

      - Change patch type:
        \verbatim
          foamDictionary constant/polyMesh/boundary \
            -entry entry0.fixedWalls.type -set patch
        \endverbatim
        This uses special parsing of Lists which stores these in the
        dictionary with keyword 'entryDDD' where DDD is the position
        in the dictionary (after ignoring the FoamFile entry).
      
      - Set nonUniform z0 to ground patch:
        \verbatim
	   setZ0 0/nut -entry boundaryField.terrain.z0 -setZ0Ground terrain
	\endverbatim

      - Set nonUniform z0 to inlet patch:
        \verbatim
	   setZ0 0/epsilon -entry boundaryField.inlet.z0 -setZ0Inlet inlet
	\endverbatim

\*---------------------------------------------------------------------------*/

//#include "fvPatchFieldMapper.H"
//#include "surfaceFields.H"
//#include <iostream>
//#include "volFields.H"
//#include "fvMesh.H"
//#include "fvBoundaryMesh.H"
//#include "polyMesh.H"
//#include "IOobjectList.H"
//#include "argList.H"
#include "Time.H"
#include "cpuTime.H"
#include "fvCFD.H"
#include "argList.H"
#include "IOobject.H"
#include "Pair.H"
#include "IFstream.H"
#include "OFstream.H"
#include "IStringStream.H"
#include "includeEntry.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//Export patch face centers' coordinates to file
void exportCf(const vectorField& face_centers, const fvBoundaryMesh& bMesh, label& patch)
{
        Info << "Writing coordinates to xyz file..." << endl;
        OFstream outputfile("Cf_coords.xyz");
        outputfile << "xyz\n";
	const labelUList inCells = bMesh[patch].faceCells();
        //Iterate through the ground faces and store to array the ground face centers coordinate
	forAll(face_centers, facei)
        {
		 scalar x_coord = face_centers[facei].x();
		 scalar y_coord = face_centers[facei].y();
		 scalar z_coord = face_centers[facei].z();
		 label bCell = inCells[facei];
		 label gfacei = bMesh[patch].start() + facei;
		 outputfile << x_coord << " " << y_coord << " " << z_coord << " " << gfacei << " " << bCell << endl;  
	}
}

//Get nocomment line. Same as Foam::getLineNoComment
string getNoCommentLine(IFstream& is)
{
       /* Same as function Foam::getLineNoComment, line 41 
	* of file in link: https://cpp.openfoam.org/v7/surfaceFormatsCore_8C_source.html#l00041*/

       string line;
       do
       {
	     is.getLine(line);
       }
       while ((line.empty() || line[0] == '#') && is.good());

       return line;
}

//Check that mtl landcover names match the ones in the app's dictionary file
void checkMtl(const fileName& MtlfileName, dictionary& landcoverDict)
{
	/* Inspired by readOBJ.C link:
	 * https://cpp.openfoam.org/v7/readOBJ_8C_source.html*/
        
	Info << "Checking landcover information in " << MtlfileName << nl << endl;

	IFstream Mtlfile(MtlfileName);

	if (!Mtlfile.good())
	{
		FatalError
		    << "Cannot read file " << MtlfileName << nl
		    << exit(FatalError);
	}
        
	int idx = 0;
	wordList landcovers;
        wordList landcoverNames(landcoverDict.toc());

	while (Mtlfile.good())
	{
		string line = getNoCommentLine(Mtlfile);
		label sz = line.size();

		if (sz && line[sz-1] == '\\') //check for block comment
		{
		     line.substr(0, sz-1);																					
		     line += getNoCommentLine(Mtlfile);
		}

	        // Read first word
		IStringStream lineStream(line);
		word fw;
		lineStream >> fw;
		if (fw == "newmtl")
		{     
		      word landcover;
		      lineStream >> landcover;
		      landcovers.append(landcover);
                      
		      //Check whether the specified landocovers is the .mtl file are found in myDict
		      auto foundc = std::find(landcoverNames.begin(), landcoverNames.end(), landcovers[idx]);
		      if (foundc == landcoverNames.end())
		      {
			      FatalError
				  << "Landcover with name " << landcovers[idx] << nl
				  << "found in the provided mtl " << MtlfileName << nl
				  << " does not match the entries in setZ0Dict. " 
				  << exit(FatalError);
		      }
		      idx++;
		}
		if (Mtlfile.eof())
		{
		      //Info << "Reached end of " << MtlfileName << endl;
		}
         }

	 if ((landcovers.empty() == true) || (landcovers.size() != landcoverNames.size()))
	 {
		 FatalError
		     << "Check that landcovers in the mtl file " << nl
		     << MtlfileName << " match the input landcovers in dictionary." << nl
		     << exit(FatalError);
	 } 
	 Info << MtlfileName << " OK. " << nl << endl;
	      //<< "Landcover information of mtl file match setZ0Dict."<< endl; 
}

//Check obj file for mtl information. Assign z0 to traingles. List corresponds to list of faces in the triSurface.
//(Potential problem if number of triangles in triSurface is different to the ones in obj...)
void z02triangle(const fileName& OBJfileName, const fileName& MtlfileName, dictionary& landcoverDict, List<double>& triZ0)
{
	/* Inspired by readObj.C Link:
	 * https://cpp.openfoam.org/v7/readOBJ_8C_source.html*/

	Info << "Checking landcover information in " << OBJfileName << nl << endl;

	IFstream OBJfile(OBJfileName);

	if (!OBJfile.good())
	{
		FatalErrorInFunction
		   << "Cannot read file " << OBJfileName
		   << exit(FatalError);
	}
        
	wordList landcoverNames(landcoverDict.toc());
        
        word landcover;	
        int counter = 1;
	int usemtl_count = 0;
	int face_count = 0;
	double val;

	while (OBJfile.good())
	{
	      string line = getNoCommentLine(OBJfile);
	      label sz = line.size();

	      if (sz && line[sz-1] == '\\')
	      {	
		    line.substr(0, sz-1);
		    line += getNoCommentLine(OBJfile);
	      }

	      // Read first word
	      IStringStream lineStream(line);
	      word fw;      
	      lineStream >> fw;

	      if (counter == 1)
              {	    
		    //Search for 'mtllib'
		    string::size_type startPos = line.find("mtllib");

                    if (startPos != string::npos)
	            {  
			//Check for space after word 'mtllib'
			startPos = line.find_first_not_of("mtllib");
                
			if (startPos == string::npos)
		        {
			     FatalError
				 << "The mtl file pathname is not specified in the " << nl
				 << "provided .obj file." << nl
				 << exit(FatalError);
			}
			
			//Get mtl name from obj file
			string mtlpath = line.substr(startPos, line.size() - startPos);
			startPos = mtlpath.rfind('/') + 1;
			string mtlname = mtlpath.substr(startPos, mtlpath.size() - startPos);

			//Get mtl name form Dict
			startPos = MtlfileName.rfind('/') + 1;
			string mtlnameDict = 
				MtlfileName.substr(startPos, MtlfileName.size() - startPos);

			//Check if mtl name in obj is same with mtl name in Dict
			string::size_type isPos = mtlnameDict.find_first_not_of(mtlname);

			if (isPos != string::npos)
			{
		              FatalError
				  << "Something is wrong in the specified .obj file." << nl
				  << "Check line " << counter << nl
				  << exit(FatalError);
		        }
                    }
		    else
		    {
			FatalError
			    << "Something is wrong in the specified .obj file." << nl
			    << "Check no comment line " << counter << nl
			    << exit(FatalError);	
                    }		    
	       }
	       else if (fw == "usemtl")						  
	       {
		    usemtl_count++;		    
		    lineStream >> landcover;
		    
		    //Check if landcover is in myDict.subdict z0_values
		    auto found_match = std::find(landcoverNames.begin(), landcoverNames.end(), landcover);
		    if (found_match == landcoverNames.end())
		    {
			    FatalError
			        << "The specified landcover in line " << counter
			       	<< " of the .obj file can not be found in the setZ0Dict. " << nl
				<< exit(FatalError); 
                    }
               }
	       else if (fw == "f")
	       { 
		       if (face_count < triZ0.size() && usemtl_count > 0)
		       {    
			    landcoverDict.lookup(landcover) >> val;
			    triZ0[face_count]  = val;
		       }
		       face_count++;
	       }
	       if (OBJfile.eof())
	       {
		       //Info << "Reached end of " << OBJfileName << endl;
	       }
	       counter++;
         }
         
	 if (usemtl_count == 0)
         {
		 FatalErrorInFunction
		     << "No entries of <usemtl landcoverName> were found." << nl
		     << exit(FatalError);
         }
	 Info << OBJfileName << " OK. " << nl << endl;
}

//Export z0 values to .txt file
void exportZ0(List<double>& z0List)
{
	OFstream outputfile("z0List.txt");
	outputfile << z0List.size() << nl << "(" << nl;
	
	forAll(z0List, index)
	{
	      outputfile << z0List[index] << nl;
	}
	outputfile << ");";

}

//Assign z0 values to face centers based on built-in octree search of nearest point
int assignZ0(triSurface& surf, const vectorField& faceCenters, List<double>& triZ0, List<double>& z0List, word& patchName)
{       
	/* Following block of code from surfaceToPoint.C
	 * Link: https://github.com/OpenFOAM/OpenFOAM-7/blob/master/src/meshTools/sets/pointSources/surfaceToPoint/surfaceToPoint.C */
        
	Info << "Assigning z0 to patch " << patchName << " face centers..." << nl << endl;

        //Construct search engine on surface
	triSurfaceSearch querySurf(surf);

	float dist = 0.8; //Needs tuning	
	int nNotAssigned = 0;
	//Box dimensions to search in octree.
	const vector span(dist, dist, dist);
    
	forAll(faceCenters, facei)
	{
		const point& pt = faceCenters[facei];
		pointIndexHit inter = querySurf.nearest(pt,span);

		if (inter.hit() && (mag(inter.hitPoint() - pt) < dist))
		{
			int triIndex = inter.index();
			z0List[facei] = triZ0[triIndex];
		}
		else
		{
			nNotAssigned += 1;
		}
	}
    
	return nNotAssigned;
}

// Enable visualisation of the result by exporting result to vtk (grid and z0 values)
void export2Vtk(const fvMesh& mesh, const word& patchName, const label& patchId, const List<double>& z0List)
{	
	Info << "Writing values to vtk file..." << endl;

	const fileName& outputDir = "constant";
	OFstream os(outputDir/"z0" + '_' + patchName + ".vtk");

	//header
	os << "# vtk DataFile Version 2.0" << nl
	   << "sampleSurface" << nl
   	   << "ASCII" << nl
	   << "DATASET POLYDATA" << nl;

	//Write vertex coords
	const pointField& points = mesh.points(); //Node coords
	const labelList& patch_pts = mesh.boundaryMesh()[patchId].meshPoints(); //Node global index 

	os << "POINTS " << patch_pts.size()
	                << " double" << nl;

	forAll(patch_pts, pI)
	{
		const point& pt = points[patch_pts[pI]];
	        os << float(pt.x()) << ' '
		   << float(pt.y()) << ' '
	   	   << float(pt.z()) << nl;
        }
	os << nl;

	//Write faces
	label nNodes = 0;
	const faceList& faces = mesh.faces(); //Face to node index

	forAll(mesh.boundary()[patchId], facei)
	{
		const label& fIdx = mesh.boundary()[patchId].start() + facei; // Face global index
		nNodes += faces[fIdx].size();
	}
	os << "POLYGONS " << mesh.boundary()[patchId].size() << ' '
	   << mesh.boundary()[patchId].size() + nNodes << nl;

	forAll(mesh.boundary()[patchId], facei)
	{
		const label& f_idx = mesh.boundary()[patchId].start() + facei;
		const face& f = faces[f_idx];
		os << f.size();

		forAll(f,fn)
		{
			auto found_pIdx = std::find(patch_pts.begin(),patch_pts.end(),f[fn]);

			if (found_pIdx != patch_pts.end())
			{
				auto pIdxInList = std::distance(patch_pts.begin(), found_pIdx);
				os << ' ' << pIdxInList;
			}
		}
		os << nl;
	}	
	
	// start writing data
	os << "CELL_DATA ";

	os << z0List.size() << nl
	   << "FIELD attributes 1" << nl
	   << "z0" << " ";

        //Write data
	os << "1 " << z0List.size() << " float" << nl;

	forAll(z0List, i)
	{
		if (i)
		{
			if (i % 10)
			{
				os << ' ';
			}
			else
			{
				os << nl;
			}
		}
		os << float(z0List[i]);
	}
	os << nl;
}
	

//Assign z0 to column based on z0 of bottom face
void assignColZ0(label faceI, const label bCell, const fvMesh& mesh, List<double>& z0, List<scalar>& heights, const labelUList& inCells)
{
	const labelList celliCells = mesh.cellCells(bCell);
	forAll(celliCells, celli)
        {       
		auto foundCelli = std::find(inCells.begin(),inCells.end(),celliCells[celli]);
		
		if (foundCelli != inCells.end())
	        {
			auto nextFaceI = std::distance(inCells.begin(), foundCelli);
			//Compare heights and skip cells when z0 is already assigned    
			if (heights[nextFaceI] > heights[faceI] && z0[nextFaceI] == 0)
			{
				z0[nextFaceI] = z0[faceI];
			 	assignColZ0(nextFaceI, inCells[nextFaceI], mesh, z0, heights, inCells);
			}
	        }
	}

}

void check_patch(const fvBoundaryMesh& bMesh, const label& patchIndex, const bool isGround)
{
	//The check does not filter out the side patch sharing cell0 index
	const label bCell0 = bMesh[patchIndex].faceCells()[0];
        //Check if patch is a ground patch or inlet patch
        if (bMesh[patchIndex].type() != "wall" && isGround == true)
	{
	        FatalErrorInFunction
	            << "The specified patch '" << bMesh[patchIndex].name() <<"' is of type '"
		    << bMesh[patchIndex].type() << "'. " << nl
		    << "A patch of type 'wall' is required." << nl		
		    << exit(FatalError);
	}
	else if (bMesh[patchIndex].type() == "wall" && isGround == false)
	{
		FatalErrorInFunction
		    << "The specified patch '" << bMesh[patchIndex].name() <<"' is of type '"
		    << bMesh[patchIndex].type() << "'. " << nl
		    << "A patch of type 'patch' is required." << nl
		    << exit(FatalError);
	}
	else if (bMesh[patchIndex].type() != "wall" && isGround == false && bCell0 != 0)
	{
		FatalErrorInFunction
		    << "The specified patch is not the inlet." << nl
		    << exit(FatalError);
	}
}


//- Read dictionary from file and return
//  Sets steam to binary mode if specified in the optional header
IOstream::streamFormat readDict(dictionary& dict, const fileName& dictFileName)
{
    IOstream::streamFormat dictFormat = IOstream::ASCII;
 
    IFstream dictFile(dictFileName);
    if (!dictFile().good())
    {   
        FatalErrorInFunction
            << "Cannot open file " << dictFileName
            << exit(FatalError, 1);
    }

    // Read the first entry from the dictionary
    autoPtr<entry> firstEntry(entry::New(dictFile()));

    // If the first entry is the "FoamFile" header dictionary
    // read and set the stream format
    if (firstEntry->isDict() && firstEntry->keyword() == "FoamFile")
    {
        dictFormat = IOstream::formatEnum(firstEntry->dict().lookup("format"));
        dictFile().format(dictFormat);
    }

    // Add the first entry to the dictionary
    dict.add(firstEntry);

    // Read and add the rest of the dictionary entries
    // preserving the "FoamFile" header dictionary if present
    dict.read(dictFile(), true);

    return dictFormat;
}


//- Converts old scope syntax to new syntax
word scope(const fileName& entryName)
{ 
    if (entryName.find(':') != string::npos)
    {
        wordList entryNames(entryName.components(':'));

        word entry(entryNames[0]);
        for (label i = 1; i < entryNames.size(); i++)
        {
            entry += word('.') + entryNames[i];
        }
        return entry;
    }
    else
    {
        return entryName;
    }
}

//- Extracts dict name and keyword
Pair<word> dictAndKeyword(const word& scopedName)
{
    string::size_type i = scopedName.find_last_of(".");
    if (i != string::npos)
    {
        return Pair<word>
        (
            scopedName.substr(0, i),
            scopedName.substr(i+1, string::npos)
        );
    }
    else
    {
        return Pair<word>("", scopedName);
    }
}

const dictionary& lookupScopedDict
(
    const dictionary& dict,
    const word& subDictName
)
{
    if (subDictName == "")
    {
        return dict;
    }
    else
    {   
        const entry* entPtr = dict.lookupScopedEntryPtr
        (
            subDictName,
            false,
            false
        );
        if (!entPtr || !entPtr->isDict())
        {
            FatalIOErrorInFunction(dict)
                << "keyword " << subDictName
                << " is undefined in dictionary "
                << dict.name() << " or is not a dictionary"
                << endl
                << "Valid keywords are " << dict.keys()
                << exit(FatalIOError);
        }
        return entPtr->dict();
    }
}


void remove(dictionary& dict, const dictionary& removeDict)
{
    forAllConstIter(dictionary, removeDict, iter)
    {
        const entry* entPtr = dict.lookupEntryPtr
        (
            iter().keyword(),
            false,
            false
        );

        if (entPtr)
        {
            if (entPtr->isDict())
            {
                if (iter().isDict())
                {
                    remove
                    (
                        const_cast<dictionary&>(entPtr->dict()),
                        iter().dict()
                    );

                    // Check if dictionary is empty
                    if (!entPtr->dict().size())
                    {
                        dict.remove(iter().keyword());
                    }
                }
            }
            else if (!iter().isDict())
            {
                if (*entPtr == iter())
                {
                    dict.remove(iter().keyword());
                }
            }
        }
    }
}

int main(int argc, char *argv[])
{   
    #include "removeCaseOptions.H"

    writeInfoHeader = false;

    argList::addNote("manipulates dictionaries");
    argList::validArgs.append("dictionary file");
    argList::addBoolOption("keywords", "list keywords");
    argList::addOption("entry", "name", "report/select the named entry");
    argList::addBoolOption
    (
        "value",
        "Print entry value"
    );
    argList::addOption
    (
        "set",
        "value",
        "Set entry value or add new entry"
    );
    argList::addOption
    (
        "add",
        "value",
        "Add a new entry"
    );
   
    argList::addBoolOption
    (
        "remove",
        "Remove the entry."
    );
   
    argList::addBoolOption
    (
        "includes",
        "List the #include/#includeIfPresent files to standard output"
    );
    argList::addBoolOption
    (
        "expand",
        "Read the specified dictionary file, expand the macros etc. and write "
        "the resulting dictionary to standard output"
    );
    argList::addBoolOption
    (
        "disableFunctionEntries",
        "Disable expansion of dictionary directives - #include, #codeStream etc"
    );
    argList::addBoolOption
    ( 
        "writeCoords",
	"Write face centers coordinates to a .xyz file"
    );
    argList::addBoolOption
    (
        "writeZ0",
        "Write z0 values to a .txt file"
    );
    argList::addOption
    (
        "setZ0Ground",
	"patchName",
	"Set or add nonUniform z0 to a ground patch. "
        "e.g 0/nut -entry boundaryField.terrain.z0 -setZ0Ground terrain. "
    );
    argList::addOption
    (
	"setZ0Inlet",
	"patchName",
	"Set or add nonUniform z0 to inlet patch. "
	"e.g 0/epsilon -entry boundaryField.inlet.z0 -setZ0Inlet inlet."
    );
    
    #include "setRootCase.H"
    #include "createTime.H"
    const bool listIncludes = args.optionFound("includes");

    if (listIncludes)
    {
        Foam::functionEntries::includeEntry::log = true;
    }

    if (args.optionFound("disableFunctionEntries"))
    {
        entry::disableFunctionEntries = true;
    }


    const fileName dictFileName(args[1]);
    dictionary dict;
    IOstream::streamFormat dictFormat = readDict(dict, dictFileName);

    bool changed = false;

    if (listIncludes)
    {
        return 0;
    }
    else if (args.optionFound("expand"))
    {
        IOobject::writeBanner(Info)
            <<"//\n// " << dictFileName << "\n//\n";
        dict.write(Info, false);
        IOobject::writeDivider(Info);

        return 0;
    }


    word entryName;
    if (args.optionReadIfPresent("entry", entryName))
    {
        word scopedName(scope(entryName));

        string newValue;
	word patch_name;
        if
        (
            args.optionReadIfPresent("set", newValue)
	 || args.optionReadIfPresent("setZ0Ground", patch_name)
	 || args.optionReadIfPresent("setZ0Inlet", patch_name)
         || args.optionReadIfPresent("add", newValue)
        )
        {
            const bool overwrite = args.optionFound("set");
	    const bool setGroundZ0 = args.optionFound("setZ0Ground");
	    const bool setInletZ0 = args.optionFound("setZ0Inlet");
            const bool merge = false; // used for add
             
            Pair<word> dAk(dictAndKeyword(scopedName));
            const dictionary& d(lookupScopedDict(dict, dAk.first()));
            
            entry* ePtr = nullptr;

	    if (setGroundZ0 || setInletZ0)
            { 
	       
	          //Access case's fvMesh object
	          //const fvMesh& mesh = runTime.lookupObject<fvMesh>("mesh"); //not possible not registered
	          ////time_db().parent.lookupObject<defaultRegion>("fvmesh");
	          //Info << mesh1.meshDir()<<endl;
		  cpuTime timer;
	          Info << "Creating instance of fvMesh..." << endl;
                  
	          fvMesh mesh
	          (
		        IOobject
		        (
			      fvMesh::defaultRegion,
			      runTime.timeName(),   //if snappy timeName = 2
			      runTime,
		              IOobject::MUST_READ,
			      IOobject::NO_WRITE,
			      false			
		        )  
	          );

		  Info << "Finished creating fvMesh instance in " << timer.cpuTimeIncrement() << " s" << endl;

	          //Check if the input patch exists in boundary mesh and return index
	          const fvBoundaryMesh& boundaryMesh = mesh.boundary();
	          label patchi = boundaryMesh.findPatchID(patch_name);
		  std::ostringstream is;
	          string s;
                  
		  Info << "Checking specified patch..." << endl;
	          if (patchi < 0)
	          {
		         FatalError
		             << "Unable to find patch " << patch_name << nl
		             << exit(FatalError);
	          }
 	          else
	          {
		        //Store all face centers of patch in vector field
		        const vectorField& Cf = boundaryMesh[patchi].Cf();
                        
			check_patch(boundaryMesh, patchi, setGroundZ0);

			IOdictionary myDict
		        (
                             IOobject
			     (
			           "setZ0Dict",
				   runTime.constant(),
				   runTime,
				   IOobject::MUST_READ
			     )
                        );

			if (args.optionFound("writeCoords"))
		        {
		             exportCf(Cf,mesh,patchi);
			}

		        //Load the input .mtl file
		        const fileName& MtlFileName(fileName(myDict.lookup("inputMtl")).expand());
                      		
			if (MtlFileName.ext() != "mtl")
			{
				FatalError
				    << "The input semantics file " << MtlFileName << nl
				    << "is not an '.mtl' file." << nl
				    << exit(FatalError);
			}
		       
		        //Check mtl file is good
		        checkMtl(MtlFileName, myDict.subDict("z0_values"));

		        //Load the input triangulation
		        const fileName& surfName_(fileName(myDict.lookup("inputFile")).expand());
                     
		        //Check if the input triangulation is in obj format
		        if (surfName_.ext() != "obj")
		        {
			       FatalError
			           << "The input geometry  " << surfName_ << nl
			           << "is in " << surfName_.ext() << " ." << nl 
			           << "An .obj file is required. " << nl
			           << exit(FatalError);
		        }

		        triSurface surf(surfName_);
		        const label surf_size = surf.faces().size(); //triangles number

		        List<double> z0tri(surf_size); //list of z0 per triangle
			z02triangle(surfName_, MtlFileName, myDict.subDict("z0_values"), z0tri);

			bool isUniform;
		        int notAssigned = 0;	
			List<double> z0patch(Cf.size(),double(0)); //list to be filled with z0 per patch face

			if (setGroundZ0)
			{
			       isUniform = false; // 0, this option is only meant for nonUniform patches else, [option] -set
			       notAssigned = assignZ0(surf, Cf, z0tri, z0patch, patch_name);
			       
			       Info << "For " << notAssigned << " face centers, out of " << boundaryMesh[patchi].size() << ", the z0 could not be assigned. " << endl; 

			       if (notAssigned != 0) //assuming that the code works perfectly
			       {
                                      FatalError
					  << "The specified patch " << patch_name << nl
					  << " does not match the " << nl
					  << "the input geometry." << nl
					  << exit(FatalError);
                               }
			       
			       //export result to vtk file
			       export2Vtk(mesh, patch_name, patchi, z0patch);

		                //From https://cpp.openfoam.org/v6/fvMesh_8C_source.html#l00246
		                /*void Foam::fvMesh::clearOut() 
			        clearAddressing();    
			        triSurface::~triSurface(); or triSurface::~clearOut();*/
                        }
                        else 
	                {                        
			       isUniform = true; // 1 
                               const labelUList& inletCells = boundaryMesh[patchi].faceCells();
			       
			       //Generate list with heights
			       List<scalar> faceHeights(Cf.size());
			       
			       forAll(faceHeights, facei)
			       {
				      faceHeights[facei] = Cf[facei].z();
			       }
			       
			       wordList checkList(Cf.size());
                              
			       forAll(boundaryMesh, patchg)
			       {                    
			              if (boundaryMesh[patchg].type() == "wall")
			              { 
					     word patchg_name = boundaryMesh[patchg].name();
			                     const labelUList& patchCells = boundaryMesh[patchg].faceCells();

					     std::stringstream z0stream;
					     double dval;
					     bool isEmpty = true; //for checking whether there is an input z0 for other patches
					     const vectorField& gCf = boundaryMesh[patchg].Cf();
					     
					     //Create list for ground patch
					     List<double> groundZ0(gCf.size());
					     
					     dictionary dnut;
					     readDict(dnut, "0/nut");
					     string patchgPath = "boundaryField." + patchg_name;
                                             const dictionary& patchgDict(lookupScopedDict(dnut, word(patchgPath)));
					     const entry* z0EntPtr = patchgDict.lookupScopedEntryPtr
					     (
					           "z0",
						   false,
						   true
					     );

					     word patchType;
					     patchgDict.lookup("type") >> patchType;
					     Info << "Type for patch in 0/nut: " << patchType << endl;

					     if (z0EntPtr) // z0 entry exists in wall patch
					     { 
					             notAssigned = assignZ0(surf, gCf, z0tri, groundZ0, patchg_name);
						     if (notAssigned != 0) //not geometry based patch
						     {
						             if (z0EntPtr->isStream())
					                     {
							             const tokenList& toks =  z0EntPtr->stream();
							             forAll(toks,t)
							             {     
									   isEmpty = true;
									   if (toks[t].isNumber())
								           {
								                  isEmpty = false;
										  z0stream << toks[t].number();
								                  z0stream >> dval;
					                                          Info << "dval " << dval << nl << endl;
                                                                           }
							             }
						              }
						     }
					     }
					     else if (!z0EntPtr && patchType != "nutkWallFunction") //excluding non ground patches
					     {       
						     Info << patchType << endl;
						     FatalIOErrorInFunction(dnut)
							   << "Can not find entry " << "z0 for wallPatch " << patchg_name
							   << exit(FatalIOError);
					     }

			                     forAll(inletCells, facei)
			                     {     
						     //bCell is global index of cell
			                             const label& bCell = inletCells[facei];
				                     auto found = std::find(patchCells.begin(), patchCells.end(), bCell);
						     if (found != patchCells.end())
				                     {
							      auto gfacei = std::distance(patchCells.begin(), found); 
							      if (notAssigned == 0)//Can cause misassignment
							      {
					                               z0patch[facei] = groundZ0[gfacei];
								       checkList[facei] = patchg_name;
						              }
							      else
						              {       
								       if (isEmpty == true)
								       {
									       FatalError
										   << "A z0 value for patch " << patchg_name << nl
										   << "has not been specified in 0/nut dictionary." << nl
										   << exit(FatalError);
								       }
								       else
								       {
								               z0patch[facei] = dval;
								       }
						              }
							      //Assign to z0 per column based on z0 of bottom face
							      assignColZ0(facei, bCell, mesh, z0patch, faceHeights, inletCells);
						      }
				               }
                                       }					   
			        }
		  		double cur;
				double old = z0patch[0];
				forAll(z0patch, facei)
			        {     
			              cur = z0patch[facei];
				      if ((cur != old) && (cur != double(0)))
                                      {
			                      isUniform = false;
			                      break;
		                      }
                                      old = cur;
			        }
		        }
		 
	                //Write z0 list to the specified file as a nonuniform list
		        if (isUniform) //setZ0Ground = 0, setZ0Inlet = (0 || 1)
	                {      
		                is << "uniform" << " " << z0patch[0] << ';';
		        }
		        else
	                {
		                is << "nonuniform List<scalar>" << nl << z0patch.size() << nl
	                        << '(' << nl;
		                for(int w=0; w<z0patch.size(); w++)
		                {
		                        is << z0patch[w] << nl;
		                }
		                is << ')' << ';';
                        }

	                if (args.optionFound("writeZ0"))
		        {
		                exportZ0(z0patch);
	                }
		 }	
		 s = is.str();
		 IStringStream str(string(dAk.second()) + " " + s);
		 ePtr = entry::New(str).ptr();
	  }
          else
          {       
                  IStringStream str(string(dAk.second()) + ' ' + newValue + ';');
                  ePtr = entry::New(str).ptr();
          }
          
          if (overwrite  || setGroundZ0  || setInletZ0)
          {
                  const_cast<dictionary&>(d).set(ePtr); //this actually clears the old ptr and uses add(entryPtr,true)     
          }
          else
          {
                  const_cast<dictionary&>(d).add(ePtr, merge); //for "add" uses add(entryPtr,false), this makes sense for checking if entry exists. Called from option <add> 
          }
          changed = true;

            // Print the changed entry
            // const entry* entPtr = dict.lookupScopedEntryPtr
            // (
            //     scopedName,
            //     false,
            //     true            // Support wildcards
            // );
            // if (entPtr)
            // {
            //     Info<< *entPtr;
            // }
        }
        else if (args.optionFound("remove"))
        {
            // Extract dictionary name and keyword
            Pair<word> dAk(dictAndKeyword(scopedName));
           
            const dictionary& d(lookupScopedDict(dict, dAk.first()));
            const_cast<dictionary&>(d).remove(dAk.second());
            changed = true;
        }
        else
        {
              const entry* entPtr = dict.lookupScopedEntryPtr
              (
                 scopedName,
                 false,
                 true            // Support wildcards
              );

              if (entPtr)
              {
                 if (args.optionFound("keywords"))
                 {
                     const dictionary& dict = entPtr->dict();
                     forAllConstIter(dictionary, dict, iter)
                     {
                         Info<< iter().keyword() << endl;
                     }
                 }
                 else
                 {
                     if (args.optionFound("value"))
                     {
                         if (entPtr->isStream())
                         {
                             const tokenList& tokens = entPtr->stream();
                             forAll(tokens, i)
                             {
                                 Info<< tokens[i];
                                 if (i < tokens.size() - 1)
                                 {
                                     Info<< token::SPACE;
                                 }
                             }
                             Info<< endl;
                         } 
                         else if (entPtr->isDict())
                         {
                             Info<< entPtr->dict();
                         }
                     }
                     else
                     {
                         Info<< *entPtr;
                     }
                  }
              }
              else
              {
                 FatalIOErrorInFunction(dict)
                    << "Cannot find entry " << entryName
                    << exit(FatalIOError, 2);
              }
         } 
    }
    else if (args.optionFound("keywords"))
    {
        forAllConstIter(dictionary, dict, iter)
        {
            Info<< iter().keyword() << endl;
        }
    }
    else
    {
        dict.write(Info, false);
    }

    if (changed)
    {  
        OFstream os(dictFileName, dictFormat);
        IOobject::writeBanner(os);
        dict.write(os, false);
        IOobject::writeEndDivider(os);
    }

    return 0;
}


// ************************************************************************* //