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
    Foam::chungKineticTransport

Description
    Real Gas Transport Package Using:
        Chung et al. (1988) model for Viscosity and Thermal Conductivity,  
        Standard Kinetic Theory model for binary diffusion coefficients.        
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
    chungKineticTransportI.H
    chungKineticTransport.C

\*---------------------------------------------------------------------------*/

#ifndef chungKineticTransport_H
#define chungKineticTransport_H

#include "PtrList.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class chungKineticTransport;

template<class Thermo>
inline chungKineticTransport<Thermo> operator+
(
    const chungKineticTransport<Thermo>&,
    const chungKineticTransport<Thermo>&
);

template<class Thermo>
inline chungKineticTransport<Thermo> operator*
(
    const scalar,
    const chungKineticTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const chungKineticTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class chungKineticTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class chungKineticTransport
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

        scalar epsilonOverKb_; 
        scalar sigma_;
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

        //- Quadratic interpolation 
        inline scalar qInterpolate(scalar t, List<scalar> X, List<scalar> Y) const;

        //- Return collision integral from tablulated data
        //- (Table IV, Monchick and Mason, Journal of Chem. Phys. 35: 1676 (1961))
        inline scalar OMEGA11(scalar TR, scalar DR) const;

public:

    // Constructors

        //- Construct from components
        inline chungKineticTransport
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
            const scalar& epsilonOverKb,
            const scalar& sigma,
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
        inline chungKineticTransport(const word&, const chungKineticTransport&);

        //- Construct from dictionary
        chungKineticTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<chungKineticTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<chungKineticTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "chungKinetic<" + Thermo::typeName() + '>';
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
        //  (mixture-averaged model) with Dij based on standard kinetic theory
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
            const List<List<scalar>>& epsilonijOverKb,
            const List<List<scalar>>& deltaij,
            const List<List<scalar>>& Mij,
            const List<List<scalar>>& sigmaij
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
        inline scalar epsilonOverKb() const;
        inline scalar sigma() const;
        inline scalar polar() const;
        inline scalar Kb() const;

        inline List<scalar> Ymd() const;
        inline List<scalar> Xmd() const;
        inline List<List<scalar>> epsilonijOverKb() const;
        inline List<List<scalar>> deltaij() const;
        inline List<List<scalar>> Mij() const;
        inline List<List<scalar>> sigmaij() const;

        
     // IO 

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline void operator=(const chungKineticTransport&); 

        inline void operator+=(const chungKineticTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend chungKineticTransport operator+ <Thermo>
        (
            const chungKineticTransport&,
            const chungKineticTransport&
        );

        friend chungKineticTransport operator* <Thermo>
        (
            const scalar,
            const chungKineticTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const chungKineticTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "chungKineticTransportI.H"

#ifdef NoRepository
    #include "chungKineticTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
