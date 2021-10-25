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

\*---------------------------------------------------------------------------*/

#include "PengRobinson.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// construct from components
template<class Specie>
inline Foam::PengRobinson<Specie>::PengRobinson
(
    const Specie& sp,
    const scalar& b_,
    const scalar& coef1_,
    const scalar& coef2_,
    const scalar& coef3_
)
:
    Specie(sp),
    b_(b_),
    coef1_(coef1_),
    coef2_(coef2_),
    coef3_(coef3_)
{}


// construct from dictionary
template<class Specie>
Foam::PengRobinson<Specie>::PengRobinson
(
    const dictionary& dict
)
:
    Specie(dict)
{
    b_ = 0.07780*RR*this->Tc_/this->Pc_ ; 

    scalar a = 0.45724*sqr(RR*this->Tc_)/this->Pc_; 

    scalar S = 0.37464 + 1.54226*this->omega_ - 0.26992*sqr(this->omega_);  
    
    coef1_ = a*sqr(1.0+S); 
 
    coef2_ = 2*a*S*(1+S)/sqrt(this->Tc_);

    coef3_ = a*sqr(S)/this->Tc_;
}


// construct as named copy
template<class Specie>
inline Foam::PengRobinson<Specie>::PengRobinson
(
    const word& name,
    const PengRobinson& srk
)
:
    Specie(name, srk),
    b_(srk.b_),
    coef1_(srk.coef1_),
    coef2_(srk.coef2_),
    coef3_(srk.coef3_)
{}


// construct and return a clone
template<class Specie>
inline Foam::autoPtr<Foam::PengRobinson <Specie>>
Foam::PengRobinson<Specie>::clone() const
{
    return autoPtr<PengRobinson<Specie>>
    (
        new PengRobinson<Specie>(*this)
    );
}


// Selector from dictionary
template<class Specie>
inline Foam::autoPtr<Foam::PengRobinson<Specie>>
Foam::PengRobinson<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<PengRobinson<Specie>>
    (
        new PengRobinson<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Specie>
void Foam::PengRobinson<Specie>::write(Ostream& os) const
{
    Specie::write(os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const PengRobinson<Specie>& srk
)
{
    srk.write(os);
    return os;
}


// ************************************************************************* //
