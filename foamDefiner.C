/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Author                                                               
    Mauro Bracconi <mauro.braccon@polimi.it>                             
    Department of Energy                                                  
    Politecnico di Milano                                                
    via La Masa, 34 - 20156 - Milano, Italy   
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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
    foamDefiner

Description
    Utility to evaluate the void fraction and the specific surface area of 
    an OpenFOAM mesh
    

\*---------------------------------------------------------------------------*/

#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
#include "fvCFD.H"
#include "multivariateScheme.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	// Solver setup (folders, mesh, etc.)
	timeSelector::addOptions();
	#include "setRootCase.H"
	#include "createTime.H"
	#include "createMesh.H"

	instantList timeDirs = timeSelector::select0(runTime, args);
	
	#include "readFoamOptions.H"

	Info<< "Time = " << runTime.timeName() << nl << endl;

	// Select foam surface patch
	label patchID = mesh.boundaryMesh().findPatchID(surfacePatch_); 
	scalar patchArea = 0;
	// if we don't have such a patch, warn the user
	if (patchID == -1) 
	{
		Info << "Failure to find patch named " << surfacePatch_ << endl;
		abort();
	}
	else  // Evaluate foam surface patch area 
	{
		const scalarField& Ap = mesh.magSf().boundaryField()[patchID];
		patchArea = sum(Ap);
		reduce(patchArea, sumOp<scalar>()); // required in parallel computing
	}

	// Evaluate volume of the mesh
	const scalarField& vol = mesh.V();
	scalar meshVolume = sum(vol);  
	reduce(meshVolume, sumOp<scalar>());    // required in parallel computing
	
	// Evaluate total volume of the bouding box
	const boundBox& boundBox = mesh.bounds();
	Vector<double> min = boundBox.min();
	Vector<double> max = boundBox.max();
	scalar bbVolume = (max[0]-min[0])*(max[1]-min[1])*(max[2]-min[2]);
	
	// Evaluate void fraction
	scalar voidFrac = 0.0;
	if( fractionType_ )
	{
		voidFrac = meshVolume/bbVolume;
	} 
	else
	{
		voidFrac = 1. - meshVolume/bbVolume;
	}

	//Write on terminal
	Info << "---------------------------------------------------" << endl;
	Info << "Catalytic area 	= " <<  patchArea << "	[m2]" << endl;
	Info << "Sv            	= " <<  patchArea/bbVolume << "	[m2/m3]" << endl;
	Info << "Void fraction 	= " <<  voidFrac << "	[-]" << endl;
	Info << "---------------------------------------------------" << endl;

	//Write on file
	autoPtr<OFstream> fmassFlux_;
	if (Pstream::parRun())
	{
		if(Pstream::master())
		{
			fmassFlux_.reset( new OFstream("../foamInfo"));
			fmassFlux_() << "Catalytic area 	= " <<  patchArea << "	[m2]" << endl;
			fmassFlux_() << "Sv            	= " <<  patchArea/bbVolume << "	[m2/m3]" << endl;
			fmassFlux_() << "Void fraction 	= " <<  voidFrac << "	[-]" << endl;
		}
	}
	else 
	{
		fmassFlux_.reset( new OFstream("foamInfo"));
		fmassFlux_() << "Catalytic area 	= " <<  patchArea << "	[m2]" << endl;
		fmassFlux_() << "Sv            	= " <<  patchArea/bbVolume << "	[m2/m3]" << endl;
		fmassFlux_() << "Void fraction 	= " <<  voidFrac << "	[-]" << endl;
	}
	
	Info << "\nEnd" << endl;
	return 0;
	
}

// ************************************************************************* //
