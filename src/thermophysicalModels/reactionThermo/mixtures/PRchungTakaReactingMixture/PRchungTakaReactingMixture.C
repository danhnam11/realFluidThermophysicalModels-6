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

#include "PRchungTakaReactingMixture.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::PRchungTakaReactingMixture<ThermoType>::PRchungTakaReactingMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    speciesTable(),
    autoPtr<chemistryReader<ThermoType>>
    (
        chemistryReader<ThermoType>::New(thermoDict, *this)
    ),
    PRchungTakaMixture<ThermoType>
    (
        thermoDict,
        *this,
        autoPtr<chemistryReader<ThermoType>>::operator()().speciesThermo(),
        mesh,
        phaseName
    ),
    PtrList<Reaction<ThermoType>>
    (
        autoPtr<chemistryReader<ThermoType>>::operator()().reactions()
    ),
    speciesComposition_
    (
        autoPtr<chemistryReader<ThermoType>>::operator()().specieComposition()
    )
{
    autoPtr<chemistryReader<ThermoType>>::clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::PRchungTakaReactingMixture<ThermoType>::read(const dictionary& thermoDict)
{}


// ************************************************************************* //
