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

#include "pureMixture.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::pureMixture<ThermoType>::pureMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicMixture(thermoDict, mesh, phaseName),
    mixture_(thermoDict.subDict("mixture")),
    ListSpeciesName_(1) //
{
    // 
    speciesData_(1);
    speciesData_.set
    (
        0,
        new ThermoType(thermoDict.subDict("mixture"))
    );
    numberOfSpecies_ = speciesData_.size();
    ListSpeciesName_[0] = "mixture";

    Y_(1);
    tmp<volScalarField> tYdefault;

    forAll(Y_, i)
    {
        IOobject header
        (
            IOobject::groupName("Y", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        );

        // check if field exists and can be read
        if (header.typeHeaderOk<volScalarField>(true))
        {
            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("Y", phaseName),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }

        else
        {
            // Read Ydefault if not already read
            if (!tYdefault.valid())
            {
                word YdefaultName(IOobject::groupName("Ydefault", phaseName));

                IOobject timeIO
                (
                    YdefaultName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                IOobject constantIO
                (
                    YdefaultName,
                    mesh.time().constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                IOobject time0IO
                (
                    YdefaultName,
                    Time::timeName(0),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                if (timeIO.typeHeaderOk<volScalarField>(true))
                {
                    tYdefault = new volScalarField(timeIO, mesh);
                }
                else if (constantIO.typeHeaderOk<volScalarField>(true))
                {
                    tYdefault = new volScalarField(constantIO, mesh);
                }
                else
                {
                    tYdefault = new volScalarField(time0IO, mesh);
                }
            }

            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("Y", phaseName),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    tYdefault()
                )
            );
        }
    }
    // 
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::pureMixture<ThermoType>::read(const dictionary& thermoDict)
{
    mixture_ = ThermoType(thermoDict.subDict("mixture"));
}


// ************************************************************************* //
