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
    Foam::sutherlandKineticTransport

Description
    Transport package using Sutherland's formula.

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

    Dynamic viscosity [kg/m.s]
    \f[
        \mu = A_s \frac{\sqrt{T}}{1 + T_s / T}
    \f]

    Mass Diffusivity is based on mixture averaged model in which 
    binary diffusion coefficients are obtained from Standard Kinetic Theory
    It has been validated against NIST and real-fluid based OPPDIF. 

    Output: 
        Dynamic Viscosity [kg/m.s]
        Thermal Conductivity [W/m.K]    
        Mass Diffusivity [m2/s]    

    by:      
    Nam Danh Nguyen - PhD.st
    Advisor: Prof. Chun Sang Yoo 
    Combustion & Propulsion Lab - Dept. of Mech. Engineering
    Ulsan Institute of Science and Technology (UNIST) - Korea

SourceFiles
    sutherlandKineticTransportI.H
    sutherlandKineticTransport.C

\*---------------------------------------------------------------------------*/

#ifndef sutherlandKineticTransport_H
#define sutherlandKineticTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class sutherlandKineticTransport;

template<class Thermo>
inline sutherlandKineticTransport<Thermo> operator+
(
    const sutherlandKineticTransport<Thermo>&,
    const sutherlandKineticTransport<Thermo>&
);

template<class Thermo>
inline sutherlandKineticTransport<Thermo> operator*
(
    const scalar,
    const sutherlandKineticTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const sutherlandKineticTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class sutherlandKineticTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class sutherlandKineticTransport
:
    public Thermo
{
    // Private data

        // Sutherland's coefficients
        scalar As_, Ts_;


    // Private Member Functions

        //- Calculate the Sutherland coefficients
        //  given two viscosities and temperatures
        inline void calcCoeffs
        (
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2
        );

        //- Read coefficient from dictionary
        scalar readCoeff(const word& coeffName, const dictionary& dict);

protected:

     //- Members for diffusivity evaluation
        scalar epsilonOverKb_;
        scalar sigma_;
        scalar miui_;
        scalar polar_;
        scalar Kb_;

        //- list of mass fraction 
        List<scalar> Ymd_;

        //- list of mole fraction 
        List<scalar> Xmd_;

        //- list of epsilonijOverKb[i][j] 
        List<List<scalar>> epsilonijOverKb_;

        //- list of deltaij[i][j] 
        List<List<scalar>> deltaij_;

        //- list of Mij[i][j] 
        List<List<scalar>> Mij_;

        //- list of sigmaij[i][j] 
        List<List<scalar>> sigmaij_;

      //- protected functions
        //- Quadratic interpolation 
        inline scalar qInterpolate(scalar t, List<scalar> X, List<scalar> Y) const;

        // Return collision integral from tablulated data
        //- (Table IV, Monchick and Mason, Journal of Chem. Phys. 35: 1676 (1961))
        inline scalar OMEGA11(scalar TR, scalar DR) const;

public:

    // Constructors

        //- Construct from components
        inline sutherlandKineticTransport
        (
            const Thermo& t,
            const scalar As,
            const scalar Ts,
            // for diffusivity
            const scalar& epsilonOverKb,
            const scalar& sigma,
            const scalar& miui,
            const scalar& polar,
            const scalar& Kb,

            const List<scalar>& Ymd,
            const List<scalar>& Xmd,
            const List<List<scalar>>& epsilonijOverKb,
            const List<List<scalar>>& deltaij,
            const List<List<scalar>>& Mij,
            const List<List<scalar>>& sigmaij
        );

        //- Construct as named copy
        inline sutherlandKineticTransport(const word&, const sutherlandKineticTransport&);

        //- Construct from dictionary
        sutherlandKineticTransport(const dictionary& dict);

        //- Construct from base thermo and dictionary
       // sutherlandKineticTransport(const Thermo& t,const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<sutherlandKineticTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<sutherlandKineticTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "sutherlandKinetic<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        // inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;

        // 
        // Species diffusivity using Kinetic model for binary diffusion coefficients
        inline scalar Dimix(label speciei, const scalar p, const scalar T) const;

        //- Update coefficients for mixture transport properties
        inline void updateTRANS
        (
            const List<scalar>& Ymd,
            const List<scalar>& Xmd,
            const List<List<scalar>>& epsilonijOverKb,
            const List<List<scalar>>& deltaij,
            const List<List<scalar>>& Mij,
            const List<List<scalar>>& sigmaij
        );

        //- Functions return member data for testing mixture  
        // for diffusivity
        inline scalar epsilonOverKb() const;
        inline scalar sigma() const;
        inline scalar miui() const;
        inline scalar polar() const;
        inline scalar Kb() const;

        inline List<scalar> Ymd() const;
        inline List<scalar> Xmd() const;
        inline List<List<scalar>> epsilonijOverKb() const;
        inline List<List<scalar>> deltaij() const;
        inline List<List<scalar>> Mij() const;
        inline List<List<scalar>> sigmaij() const;

        //

    // Member operators

        inline void operator=(const sutherlandKineticTransport&);

        inline void operator+=(const sutherlandKineticTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend sutherlandKineticTransport operator+ <Thermo>
        (
            const sutherlandKineticTransport&,
            const sutherlandKineticTransport&
        );

        friend sutherlandKineticTransport operator* <Thermo>
        (
            const scalar,
            const sutherlandKineticTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const sutherlandKineticTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sutherlandKineticTransportI.H"

#ifdef NoRepository
    #include "sutherlandKineticTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
