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
    Foam::chungTakaTransport

Description
    Real Gas Transport Package Using:
        Chung et al. (1988) model for Viscosity and Thermal Conductivity,  
        Fuller (1987) and Takahashi et al. (1975) model for binary diffusion coefficients. 
        Mass Diffusivity is based on mixture averaged model.
    ----------------------------
    has been validated against NIST and real-fluid based OPPDIF. 

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
    chungTakaTransportI.H
    chungTakaTransport.C

\*---------------------------------------------------------------------------*/

#ifndef chungTakaTransport_H
#define chungTakaTransport_H

#include "PtrList.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class chungTakaTransport;

template<class Thermo>
inline chungTakaTransport<Thermo> operator+
(
    const chungTakaTransport<Thermo>&,
    const chungTakaTransport<Thermo>&
);

template<class Thermo>
inline chungTakaTransport<Thermo> operator*
(
    const scalar,
    const chungTakaTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const chungTakaTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class chungTakaTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class chungTakaTransport
:
    public Thermo
{
    // Private data

    // Private function

    // Protected data

protected: 
        
    // All following members need to be updated for mixture at mixture classes
     //- Members for Viscosity and Thermal Conductivity evaluation
        //- return sigma in eq. (4) - Chung's Model 
        scalar sigmaM_;   

        //- return epsilon/k in eq. (5) - Chung's Model
        scalar epsilonkM_;  

        //- molecular weight
        scalar MM_;
  
        //- critical volume 
        scalar VcM_;  

        //- critical temperature
        scalar TcM_;  

        //- omega
        scalar omegaM_;  

        //- miui
        scalar miuiM_;   

        //- kappai 
        scalar kappaiM_;  
 
     //- Members for diffusivity evaluation
        //- list of mass fraction 
        List<scalar> Ymd_; 

        //- list of mole fraction 
        List<scalar> Xmd_; 

        //- list of Tcmd[i][j] 
        List<List<scalar>> Tcmd_; 

        //- list of Pcmd[i][j] 
        List<List<scalar>> Pcmd_; 

        //- list of Mmd[i][j] 
        List<List<scalar>> Mmd_; 

        //- list of sigmd[i][j] 
        List<List<scalar>> sigmd_; 
 
     //- protected function
        inline void calculate
        (
            scalar p,
            scalar T,
            scalar& Tstar,
            scalar& mur,
            scalar& eta0,
            scalar& Y,
            scalar& G1 
        ) const;

   
public:

    // Constructors

        //- Construct from components
        inline chungTakaTransport
        (
            const Thermo& t,
            const scalar& sigmaM,
            const scalar& epsilonkM,
            const scalar& MM,
            const scalar& VcM,
            const scalar& TcM,
            const scalar& omegaM,
            const scalar& miuiM,
            const scalar& kappaiM,
            // for diffusivity
            const List<scalar>& Ymd,
            const List<scalar>& Xmd,
            const List<List<scalar>>& Tcmd,
            const List<List<scalar>>& Pcmd,
            const List<List<scalar>>& Mmd,
            const List<List<scalar>>& sigmd
        );

        //- Construct as named copy
        inline chungTakaTransport(const word&, const chungTakaTransport&);

        //- Construct from dictionary
        chungTakaTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<chungTakaTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<chungTakaTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "chungTaka<" + Thermo::typeName() + '>';
        }

     // Fundamental properties 

        //- Return Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Return Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Return Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        //- Return coefficient Phi based on Takahashi model to calculate
        //  the Binary Diffusion Coefficients in Fuller model
        inline scalar phi(scalar PRR, scalar TRR) const;

        //- Return mass Diffusivity of specie ith [m^2/s] using Fick's law
        //  (mixture-averaged model) with Dij based on Fuller and Takahashi models
        inline scalar Dimix(label speciei, const scalar p, const scalar T) const;

        //- Update coefficients for mixture transport properties
        inline void updateTRANS
        (
            const scalar& sigmaM, 
            const scalar& epsilonkM,
            const scalar& MM,
            const scalar& VcM,
            const scalar& TcM,
            const scalar& omegaM,
            const scalar& miuiM,
            const scalar& kappaiM,
            const List<scalar>& Ymd, 
            const List<scalar>& Xmd, 
            const List<List<scalar>>& Tcmd,
            const List<List<scalar>>& Pcmd,
            const List<List<scalar>>& Mmd,
            const List<List<scalar>>& sigmd
        );

        //- Functions return member data for testing mixture  
        inline scalar sigmaM() const;
        inline scalar epsilonkM() const;
        inline scalar MM() const;
        inline scalar VcM() const;
        inline scalar TcM() const;
        inline scalar omegaM() const;
        inline scalar miuiM() const;
        inline scalar kappaiM() const;
        // for diffusivity
        inline List<scalar> Ymd() const;
        inline List<scalar> Xmd() const;
        inline List<List<scalar>> Tcmd() const;
        inline List<List<scalar>> Pcmd() const;
        inline List<List<scalar>> Mmd() const;
        inline List<List<scalar>> sigmd() const;

        
     // IO 

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline void operator=(const chungTakaTransport&); 

        inline void operator+=(const chungTakaTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend chungTakaTransport operator+ <Thermo>
        (
            const chungTakaTransport&,
            const chungTakaTransport&
        );

        friend chungTakaTransport operator* <Thermo>
        (
            const scalar,
            const chungTakaTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const chungTakaTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "chungTakaTransportI.H"

#ifdef NoRepository
    #include "chungTakaTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
