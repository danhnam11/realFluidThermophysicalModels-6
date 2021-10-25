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
    Foam::PRchungTakaMixture

Description
    Foam::PRchungTakaMixture

    calculate real gas properties based on new mixture model
    that can account for various mixing rules:
       EoS          : PR 
       thermodynamic: rfJanaf
       transport    : ChungTaka

    by:      
    Nam Danh Nguyen - PhD.st
    Advisor: Prof. Chun Sang Yoo 
    Combustion & Propulsion Lab - Dept. of Mech. Engineering
    Ulsan Institute of Science and Technology (UNIST) - Korea

SourceFiles
    PRchungTakaMixture.C

\*---------------------------------------------------------------------------*/

#ifndef PRchungTakaMixture_H
#define PRchungTakaMixture_H

#include "basicSpecieMixture.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class PRchungTakaMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class PRchungTakaMixture
:
    public basicSpecieMixture
{
    // Private data
 
        //- Species data
        PtrList<ThermoType> speciesData_; 

        //- Temporary storage for the cell/face mixture thermo data
        mutable ThermoType mixture_;

        //- Temporary storage for the volume weighted
        //  cell/face mixture thermo data
        mutable ThermoType mixtureVol_;


    // Private Member Functions

        //- Construct the species data from the given dictionary and return the
        //  data for the first specie to initialise the mixture thermo data
        const ThermoType& constructSpeciesData(const dictionary& thermoDict);

        //- Correct the mass fractions to sum to 1
        void correctMassFractions();

        //- Calculate parameters/coefficients in real-fluid models
        void calculateRealGas
        (
	    List<scalar> X,
	    scalar& bM,
	    scalar& coef1,
	    scalar& coef2,
	    scalar& coef3,
	    scalar& sigmaM,
	    scalar& epsilonkM,
	    scalar& MM,
	    scalar& VcM,
	    scalar& TcM,
	    scalar& omegaM, 
	    scalar& miuiM,
	    scalar& kappaiM
        ) const;     

        //- Construct as copy (not implemented)
        PRchungTakaMixture(const PRchungTakaMixture<ThermoType>&);

protected:

        //- check real fluid mode 
//        bool realFluidMode_;
  
        //- Number of species in mixture
        label numberOfSpecies_; 

        //- Return list of species name
        List<word> ListSpeciesName_; 

     //- Members to store precalculations in real-fluid models
        //- List of species molecular weight
        List<scalar> ListW_;
        //- For PR
        List<scalar> BM_;
        List<List<scalar>> COEF1_;
        List<List<scalar>> COEF2_;
        List<List<scalar>> COEF3_; 
        //- For viscosity and thermal conductivity in Chung's model
        List<List<scalar>> SIGMA3M_;
        List<List<scalar>> EPSILONKM0_;
        List<List<scalar>> OMEGAM0_;
        List<List<scalar>> MM0_;
        List<List<scalar>> MIUIM0_;
        List<List<scalar>> KAPPAIM_;
        //- For mass diffusivity in Takahashi's correction
        List<scalar> ListTc_;
        List<scalar> ListPc_;
        List<List<scalar>> MMD_;
        List<List<scalar>> SIGMD_;

public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary, specie names, thermo database,
        //  mesh and phase name
        PRchungTakaMixture
        (
            const dictionary&,
            const wordList& specieNames,
            const HashPtrTable<ThermoType>& thermoData,
            const fvMesh&,
            const word&
        );

        //- Construct from dictionary, mesh and phase name
        PRchungTakaMixture(const dictionary&, const fvMesh&, const word&);


    //- Destructor
    virtual ~PRchungTakaMixture()
    {}


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "PRchungTakaMixture<" + ThermoType::typeName() + '>';
        }

        const ThermoType& cellMixture(const label celli) const;

        const ThermoType& patchFaceMixture
        (
            const label patchi,
            const label facei
        ) const;

        const ThermoType& cellVolMixture
        (
            const scalar p,
            const scalar T,
            const label celli
        ) const;

        const ThermoType& patchFaceVolMixture
        (
            const scalar p,
            const scalar T,
            const label patchi,
            const label facei
        ) const;

        //- Return the raw specie thermodynamic data
        const PtrList<ThermoType>& speciesData() const
        {
            return speciesData_;
        }

        //- Read dictionary
        void read(const dictionary&);

        //- Return thermo based on index
        inline const ThermoType& getLocalThermo(const label speciei) const
        {
            return speciesData_[speciei];
        }


        //- Return the list of specie name
        const List<word>& ListSpeciesName() const
        {
            return ListSpeciesName_;
        }

        //- Return number of species in the mixture
        inline const label& numberOfSpecies() const
        {
            return numberOfSpecies_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PRchungTakaMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
