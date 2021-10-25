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

#include "soaveRedlichKwong.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::soaveRedlichKwong<Specie>::soaveRedlichKwong
(
    const dictionary& dict
)
:
    Specie(dict)
{
    b_ = 0.08664*RR*this->Tc_/this->Pc_ ; 

    scalar a = 0.42747*sqr(RR*this->Tc_)/this->Pc_; 
    scalar S = 0.48508 + 1.5517*this->omega_ - 0.15613*sqr(this->omega_);  
    
    coef1_ = a*sqr(1.0+S); 

    coef2_ = (0.42747*sqr(RR*this->Tc_)/this->Pc_)*
          2*(0.48508 + 1.5517*this->omega_ - 0.15613*sqr(this->omega_))* 
          (1+(0.48508 + 1.5517*this->omega_ - 0.15613*sqr(this->omega_)))/
          sqrt(this->Tc_);  

    coef3_ = (0.42747*sqr(RR*this->Tc_)/this->Pc_)*
             sqr(0.48508 + 1.5517*this->omega_ - 0.15613*sqr(this->omega_))/this->Tc_; 
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Specie>
void Foam::soaveRedlichKwong<Specie>::write(Ostream& os) const
{
    Specie::write(os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const soaveRedlichKwong<Specie>& srk
)
{
    srk.write(os);
    return os;
}


// ************************************************************************* //
