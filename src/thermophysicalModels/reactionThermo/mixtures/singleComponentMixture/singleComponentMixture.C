/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "singleComponentMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::singleComponentMixture<ThermoType>::singleComponentMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicSpecieMixture(thermoDict, wordList(), mesh, phaseName),
    thermo_(thermoDict.subDict("mixture")),
    ListSpeciesName_(1) //
{
    // 
    speciesData_(1);
    speciesData_.set
       (
          0,
          new ThermoType("mixture", thermoDict.subDict("mixture"))
       );
    numberOfSpecies_ = speciesData_.size();
    ListSpeciesName_[0] = "mixture";
    //
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::singleComponentMixture<ThermoType>::~singleComponentMixture()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::singleComponentMixture<ThermoType>::read
(
    const dictionary& thermoDict
)
{
    thermo_ = ThermoType(thermoDict.subDict("mixture"));
}


// ************************************************************************* //
