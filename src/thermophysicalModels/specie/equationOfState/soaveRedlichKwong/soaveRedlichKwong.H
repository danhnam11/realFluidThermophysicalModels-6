/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2017 OpenFOAM Foundation
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
    Foam::soaveRedlichKwong

Description
    soaveRedlichKwong gas equation of state.
    
    have been validated against NIST and real-fluid based OPPDIF. 

    by:      
    Jae Won Shim - Ms.st, Ki Sung Jung PhD.st, and Nam Danh Nguyen - PhD.st
    Advisor: Prof. Chun Sang Yoo 
    Combustion & Propulsion Lab - Dept. of Mech. Engineering
    Ulsan Institute of Science and Technology (UNIST) - Korea

SourceFiles
    soaveRedlichKwongI.H
    soaveRedlichKwong.C

\*---------------------------------------------------------------------------*/

#ifndef soaveRedlichKwong_H
#define soaveRedlichKwong_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class soaveRedlichKwong;

template<class Specie>
inline soaveRedlichKwong<Specie> operator+
(
    const soaveRedlichKwong<Specie>&,
    const soaveRedlichKwong<Specie>&
);

template<class Specie>
inline soaveRedlichKwong<Specie> operator*
(
    const scalar,
    const soaveRedlichKwong<Specie>&
);

template<class Specie>
inline soaveRedlichKwong<Specie> operator==
(
    const soaveRedlichKwong<Specie>&,
    const soaveRedlichKwong<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const soaveRedlichKwong<Specie>&
);



/*---------------------------------------------------------------------------*\
                           Class soaveRedlichKwong Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class soaveRedlichKwong
:
    public Specie
{


protected:

        //- Parameters in SRK model
        scalar b_;
        //- New decomposed parameters in aAlpha eq.
        scalar coef1_;
        scalar coef2_;
        scalar coef3_;
        

public:

    // Constructors

        //- Construct from components
        inline soaveRedlichKwong
        (
            const Specie& sp,
            const scalar& b_,
            const scalar& coef1_,
            const scalar& coef2_,
            const scalar& coef3_
        );

        //- Construct from dictionary
        soaveRedlichKwong(const dictionary& dict);

        //- Construct as named copy
        inline soaveRedlichKwong(const word& name, const soaveRedlichKwong&);

        //- Construct and return a clone
        inline autoPtr<soaveRedlichKwong> clone() const;

        // Selector from dictionary
        inline static autoPtr<soaveRedlichKwong> New
        (
            const dictionary& dict
        );


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "soaveRedlichKwong<" + word(Specie::typeName_()) + '>';
        }

        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy departure [J/kg]
            inline scalar H(const scalar p, const scalar T) const;

            //- Return Cp departure [J/(kg K]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            inline scalar CpMCv(scalar p, scalar T) const;

            //- Update coeffcients for mixture
            void updateEoS
            ( 
                const scalar& b, 
                const scalar& coef1, 
                const scalar& coef2, 
                const scalar& coef3
            ) 
            { 
                b_ = b; 
                coef1_ = coef1; 
                coef2_ = coef2; 
                coef3_ = coef3; 
            } 

            // Functions to check parameters
            inline scalar bM() const;
            inline scalar coef1() const;
            inline scalar coef2() const;
            inline scalar coef3() const;


    // IO

            //- Write to Ostream
            void write(Ostream& os) const;

    // Member operators

        inline void operator+=(const soaveRedlichKwong&);
        inline void operator*=(const scalar);


    // Friend operators

        friend soaveRedlichKwong operator+ <Specie>
        (
            const soaveRedlichKwong&,
            const soaveRedlichKwong&
        );

        friend soaveRedlichKwong operator* <Specie>
        (
            const scalar s,
            const soaveRedlichKwong&
        );

        friend soaveRedlichKwong operator== <Specie>
        (
            const soaveRedlichKwong&,
            const soaveRedlichKwong&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const soaveRedlichKwong&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "soaveRedlichKwongI.H"

#ifdef NoRepository
    #include "soaveRedlichKwong.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
