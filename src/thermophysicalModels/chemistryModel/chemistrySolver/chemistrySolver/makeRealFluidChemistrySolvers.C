/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

    Make chemistry solvers for real-fluid thermophysicalModels library

    by:      
    Nam Danh Nguyen - PhD.st
    Advisor: Prof. Chun Sang Yoo 
    Combustion & Propulsion Lab - Dept. of Mech. Engineering
    Ulsan Institute of Science and Technology (UNIST) - Korea 	
	

\*---------------------------------------------------------------------------*/

#include "makeRealFluidChemistrySolverTypes.H"

#include "thermoPhysicsTypes.H"
#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry solvers based on sensibleEnthalpy
    makeSRKchungTakaChemistrySolverTypes(psiReactionThermo, chungTakaRealJsrkHThermoPhysics);
    makeSRKchungTakaChemistrySolverTypes(rhoReactionThermo, chungTakaRealJsrkHThermoPhysics);

    makeSRKchungKineticChemistrySolverTypes(psiReactionThermo, chungKineticRealJsrkHThermoPhysics);
    makeSRKchungKineticChemistrySolverTypes(rhoReactionThermo, chungKineticRealJsrkHThermoPhysics);

    makeidKineticChemistrySolverTypes(psiReactionThermo, idKineticHThermoPhysics);
    makeidKineticChemistrySolverTypes(rhoReactionThermo, idKineticHThermoPhysics);

    makePRchungTakaChemistrySolverTypes(psiReactionThermo, chungTakaRealJprHThermoPhysics);
    makePRchungTakaChemistrySolverTypes(rhoReactionThermo, chungTakaRealJprHThermoPhysics);

    makePRchungKineticChemistrySolverTypes(psiReactionThermo, chungKineticRealJprHThermoPhysics);
    makePRchungKineticChemistrySolverTypes(rhoReactionThermo, chungKineticRealJprHThermoPhysics);


    // Chemistry solvers based on sensibleInternalEnergy
    makeSRKchungTakaChemistrySolverTypes(psiReactionThermo, chungTakaRealJsrkEThermoPhysics);
    makeSRKchungTakaChemistrySolverTypes(rhoReactionThermo, chungTakaRealJsrkEThermoPhysics);

    makeSRKchungKineticChemistrySolverTypes(psiReactionThermo, chungKineticRealJsrkEThermoPhysics);
    makeSRKchungKineticChemistrySolverTypes(rhoReactionThermo, chungKineticRealJsrkEThermoPhysics);

    makeidKineticChemistrySolverTypes(psiReactionThermo, idKineticEThermoPhysics);
    makeidKineticChemistrySolverTypes(rhoReactionThermo, idKineticEThermoPhysics);

    makePRchungTakaChemistrySolverTypes(psiReactionThermo, chungTakaRealJprEThermoPhysics);
    makePRchungTakaChemistrySolverTypes(rhoReactionThermo, chungTakaRealJprEThermoPhysics);

    makePRchungKineticChemistrySolverTypes(psiReactionThermo, chungKineticRealJprEThermoPhysics);
    makePRchungKineticChemistrySolverTypes(rhoReactionThermo, chungKineticRealJprEThermoPhysics);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
