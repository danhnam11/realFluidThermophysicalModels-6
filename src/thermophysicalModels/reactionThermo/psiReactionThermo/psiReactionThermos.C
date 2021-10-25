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

\*---------------------------------------------------------------------------*/

#include "makeReactionThermo.H"

#include "psiReactionThermo.H"
#include "hePsiThermo.H"

#include "specie.H"
#include "rfSpecie.H" //
#include "perfectGas.H"
#include "soaveRedlichKwong.H" //
#include "PengRobinson.H" //
#include "hConstThermo.H"
#include "janafThermo.H"
#include "rfJanafThermo.H"  //
#include "sensibleEnthalpy.H"
#include "thermo.H"
#include "constTransport.H"
#include "sutherlandTransport.H"
#include "sutherlandKineticTransport.H" //
#include "chungTakaTransport.H" //
#include "chungKineticTransport.H" //

#include "homogeneousMixture.H"
#include "inhomogeneousMixture.H"
#include "veryInhomogeneousMixture.H"
#include "multiComponentMixture.H"
#include "reactingMixture.H"
#include "SRKchungTakaMixture.H" //
#include "SRKchungTakaReactingMixture.H" //
#include "SRKchungKineticMixture.H" //
#include "SRKchungKineticReactingMixture.H" //
#include "PRchungTakaMixture.H" //
#include "PRchungTakaReactingMixture.H" //
#include "PRchungKineticMixture.H" //
#include "PRchungKineticReactingMixture.H" //
#include "idKineticMixture.H" //
#include "idKineticReactingMixture.H" //
#include "singleStepReactingMixture.H"
#include "singleComponentMixture.H"

#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// constTransport, hConstThermo

makeReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    homogeneousMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    inhomogeneousMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    veryInhomogeneousMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);


// sutherlandTransport, hConstThermo

makeReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    homogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    inhomogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    veryInhomogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);


// sutherlandTransport, janafThermo

makeReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    homogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas,
    specie
);

makeReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    inhomogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas,
    specie
);

makeReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    veryInhomogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas,
    specie
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// real fluid mixture thermo for sensible enthalpy

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    idKineticMixture,
    idKineticHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    SRKchungTakaMixture,
    chungTakaRealJsrkHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    SRKchungKineticMixture,
    chungKineticRealJsrkHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    PRchungTakaMixture,
    chungTakaRealJprHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    PRchungKineticMixture,
    chungKineticRealJprHThermoPhysics
);


// real fluid mixture thermo for internal energy 
makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    idKineticMixture,
    idKineticEThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    SRKchungTakaMixture,
    chungTakaRealJsrkEThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    SRKchungKineticMixture,
    chungKineticRealJsrkEThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    PRchungTakaMixture,
    chungTakaRealJprEThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    PRchungKineticMixture,
    chungKineticRealJprEThermoPhysics
);


// real fluid mixture reaction thermo for sensible enthalpy
makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    idKineticReactingMixture,
    idKineticHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    SRKchungTakaReactingMixture,
    chungTakaRealJsrkHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    SRKchungKineticReactingMixture,
    chungKineticRealJsrkHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    PRchungTakaReactingMixture,
    chungTakaRealJprHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    PRchungKineticReactingMixture,
    chungKineticRealJprHThermoPhysics
);

// real fluid mixture reaction thermo for internal energy 
makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    idKineticReactingMixture,
    idKineticEThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    SRKchungTakaReactingMixture,
    chungTakaRealJsrkEThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    SRKchungKineticReactingMixture,
    chungKineticRealJsrkEThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    PRchungTakaReactingMixture,
    chungTakaRealJprEThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    PRchungKineticReactingMixture,
    chungKineticRealJprEThermoPhysics
);
// ------


// Multi-component thermo for sensible enthalpy

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    multiComponentMixture,
    constGasHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    multiComponentMixture,
    gasHThermoPhysics
);

// Multi-component thermo for internal energy

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    multiComponentMixture,
    constGasEThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    multiComponentMixture,
    gasEThermoPhysics
);


// Reaction thermo for sensible enthalpy

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    reactingMixture,
    constGasHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    reactingMixture,
    gasHThermoPhysics
);


// Single-step reaction thermo for sensible enthalpy

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    singleStepReactingMixture,
    gasHThermoPhysics
);


// Reaction thermo for internal energy

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    reactingMixture,
    constGasEThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    reactingMixture,
    gasEThermoPhysics
);


// Single-step reaction thermo for internal energy

makeThermoPhysicsReactionThermos
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    singleStepReactingMixture,
    gasEThermoPhysics
);


// Single-component thermo for sensible enthalpy

makeThermoPhysicsReactionThermo
(
    psiReactionThermo,
    hePsiThermo,
    singleComponentMixture,
    constGasHThermoPhysics
);

makeThermoPhysicsReactionThermo
(
    psiReactionThermo,
    hePsiThermo,
    singleComponentMixture,
    gasHThermoPhysics
);


// Single-component thermo for internal energy

makeThermoPhysicsReactionThermo
(
    psiReactionThermo,
    hePsiThermo,
    singleComponentMixture,
    constGasEThermoPhysics
);

makeThermoPhysicsReactionThermo
(
    psiReactionThermo,
    hePsiThermo,
    singleComponentMixture,
    gasEThermoPhysics
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
