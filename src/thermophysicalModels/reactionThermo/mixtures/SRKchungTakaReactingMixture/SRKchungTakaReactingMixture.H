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

Class
    Foam::SRKchungTakaReactingMixture

Description
    Foam::SRKchungTakaReactingMixture

    calculate real gas properties for reacting mixture based on new mixture model
    that can account for various mixing rules:
       EoS          : SRK 
       thermodynamic: rfJanaf
       transport    : ChungTaka

    by:      
    Nam Danh Nguyen - PhD.st
    Advisor: Prof. Chun Sang Yoo 
    Combustion & Propulsion Lab - Dept. of Mech. Engineering
    Ulsan Institute of Science and Technology (UNIST) - Korea

SourceFiles
    SRKchungTakaReactingMixture.C

\*---------------------------------------------------------------------------*/

#ifndef SRKchungTakaReactingMixture_H
#define SRKchungTakaReactingMixture_H

#include "speciesTable.H"
#include "chemistryReader.H"
#include "SRKchungTakaMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class SRKchungTakaReactingMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class SRKchungTakaReactingMixture
:
    public speciesTable,
    public autoPtr<chemistryReader<ThermoType>>,
    public SRKchungTakaMixture<ThermoType>,
    public PtrList<Reaction<ThermoType>>
{
    // Private member data

        //- Table of species composition
        speciesCompositionTable speciesComposition_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        SRKchungTakaReactingMixture(const SRKchungTakaReactingMixture&);

        //- Disallow default bitwise assignment
        void operator=(const SRKchungTakaReactingMixture&);


public:

    //- The type of thermo package this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary, mesh and phase name
        SRKchungTakaReactingMixture(const dictionary&, const fvMesh&, const word&);


    //- Destructor
    virtual ~SRKchungTakaReactingMixture()
    {}


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "SRKchungTakaReactingMixture<" + ThermoType::typeName() + '>';
        }

        //- Read dictionary
        void read(const dictionary&);

        label size() const
        {
            return PtrList<Reaction<ThermoType>>::size();
        }

        Reaction<ThermoType>& operator[](const label i)
        {
            return PtrList<Reaction<ThermoType>>::operator[](i);
        }

        const Reaction<ThermoType>& operator[](const label i) const
        {
            return PtrList<Reaction<ThermoType>>::operator[](i);
        }

        //- Table of species composition
        const speciesCompositionTable& specieComposition() const
        {
            return speciesComposition_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SRKchungTakaReactingMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
