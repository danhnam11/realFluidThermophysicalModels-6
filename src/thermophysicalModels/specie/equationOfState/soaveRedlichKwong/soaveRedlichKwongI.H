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
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::soaveRedlichKwong<Specie>::soaveRedlichKwong
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::soaveRedlichKwong<Specie>::soaveRedlichKwong
(
    const word& name,
    const soaveRedlichKwong& srk
)
:
    Specie(name, srk),
    b_(srk.b_),
    coef1_(srk.coef1_),
    coef2_(srk.coef2_),
    coef3_(srk.coef3_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::soaveRedlichKwong <Specie>>
Foam::soaveRedlichKwong<Specie>::clone() const
{
    return autoPtr<soaveRedlichKwong<Specie>>
    (
        new soaveRedlichKwong<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::soaveRedlichKwong<Specie>>
Foam::soaveRedlichKwong<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<soaveRedlichKwong<Specie>>
    (
        new soaveRedlichKwong<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Specie>
inline Foam::scalar Foam::soaveRedlichKwong<Specie>::rho
(
    scalar p,
    scalar T
) const
{
    const scalar Z = this->Z(p, T);
    return p/(Z*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::soaveRedlichKwong<Specie>::H(scalar p, scalar T) const
{
    scalar B = b_*p/(RR*T); 
    if(B <=0 ) { B = 1e-16; } 

    const scalar aAlpha = coef1_ - coef2_*sqrt(T) + coef3_*T;
    const scalar daAlpha = -coef2_/(2*sqrt(T))+coef3_;

    const scalar Z = this->Z(p, T);
    if (B == -Z){
    return (RR*T*(Z - 1))/this->W();
    }
    else {
    return (RR*T*(Z - 1)+(T*daAlpha-aAlpha)/b_ *log((Z + B)/Z))/this->W();
    }
}


template<class Specie>
inline Foam::scalar Foam::soaveRedlichKwong<Specie>::Cp(scalar p, scalar T) const
{
    const scalar aAlpha = coef1_ - coef2_*sqrt(T) + coef3_*T;
    const scalar daAlpha = -coef2_/(2*sqrt(T))+coef3_;
    const scalar ddaAlpha = coef2_/(4*T*sqrt(T));
    const scalar A = aAlpha*p/sqr(RR*T);
    scalar B = b_*p/(RR*T); 
    if(B <=0 ) { B = 1e-16; }

    const scalar Z = this->Z(p, T);

    const scalar M = (sqr(Z) + B*Z)/(Z - B);
    const scalar N = daAlpha*B/(b_*RR);
 
    return
    (
        (T/b_)*ddaAlpha*log((Z + B)/Z)
      + RR*sqr(M - N)/(sqr(M) - A*(2*Z + B))
      - RR
    )/this->W();
   

}


// This function has not been validated yet
template<class Specie>
inline Foam::scalar Foam::soaveRedlichKwong<Specie>::S
(
    scalar p,
    scalar T
) const
{
    const scalar daAlpha = -coef2_/(2*sqrt(T))+coef3_;
    scalar B = b_*p/(RR*T); 
    if(B <=0 ) { B = 1e-16; }

    const scalar Z = this->Z(p, T);

    return
        (
        RR*log(Z - B)
        -daAlpha/b_*log((Z + B)/Z)
        )/this->W();
}


template<class Specie>
inline Foam::scalar Foam::soaveRedlichKwong<Specie>::psi
(
    scalar p,
    scalar T
) const
{
    const scalar Z = this->Z(p, T);

    return 1.0/(Z*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::soaveRedlichKwong<Specie>::Z
(
    scalar p,
    scalar T
) const
{
    const scalar aAlpha = coef1_ - coef2_*sqrt(T) + coef3_*T;

    const scalar A = aAlpha*p/sqr(RR*T);
    const scalar B = b_*p/(RR*T);

    const scalar a2 = - 1.0;
    const scalar a1 = A - B - sqr(B);
    const scalar a0 = -A*B;

    const scalar Q = (3*a1 - a2*a2)/9.0;
    const scalar Rl = (9*a2*a1 - 27*a0 - 2*a2*a2*a2)/54.0;

    const scalar Q3 = Q*Q*Q;
    const scalar D = Q3 + Rl*Rl;

    scalar root = -1;

    if (D <= 0)
    {
        const scalar th = ::acos(Rl/sqrt(-Q3));
        const scalar qm = 2*sqrt(-Q);
        const scalar r1 = qm*cos(th/3.0) - a2/3.0;
        const scalar r2 =
            qm*cos((th + 2*constant::mathematical::pi)/3.0) - a2/3.0;
        const scalar r3 =
            qm*cos((th + 4*constant::mathematical::pi)/3.0) - a2/3.0;

        root = max(r1, max(r2, r3));
    }
    else
    {
        // One root is real
        const scalar D05 = sqrt(D);
        scalar S = 0;
    
        if (Rl + D05 < 0)
        {
            S = -pow(mag(Rl + D05), 1.0/3.0);
        }
        else
        {
            S = pow(Rl + D05, 1.0/3.0);
        }
        scalar Tl = 0;
        if (D05 > Rl)
        {
            Tl = -pow(mag(Rl - D05), 1.0/3.0);
        }
        else
        {
            Tl = pow(Rl - D05, 1.0/3.0);
        }

        root = S + Tl - a2/3.0;
    }
    
    return root;
}


template<class Specie>
inline Foam::scalar Foam::soaveRedlichKwong<Specie>::CpMCv
(
    scalar p,
    scalar T
) const
{
    const scalar aAlpha = coef1_ - coef2_*sqrt(T) + coef3_*T;
    const scalar daAlpha = -coef2_/(2*sqrt(T))+coef3_;

    const scalar A = aAlpha*p/sqr(RR*T);
    scalar B = b_*p/(RR*T); 
    if(B <=0 ) { B = 1e-16; } 

    const scalar Z = this->Z(p, T);

    const scalar M = (sqr(Z) + B*Z)/(Z - B);
    const scalar N = daAlpha*B/(b_*RR);

    return this->R()*sqr(M - N)/(sqr(M) - A*(2*Z + B));
}


// Functions returns member data for checking 
template<class Specie>
inline Foam::scalar Foam::soaveRedlichKwong<Specie>::bM() const
{
    return b_; 
}

template<class Specie>
inline Foam::scalar Foam::soaveRedlichKwong<Specie>::coef1() const
{
    return coef1_; 
}

template<class Specie>
inline Foam::scalar Foam::soaveRedlichKwong<Specie>::coef2() const
{
    return coef2_; 
}

template<class Specie>
inline Foam::scalar Foam::soaveRedlichKwong<Specie>::coef3() const
{
    return coef3_; 
}



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::soaveRedlichKwong<Specie>::operator+=
(
    const soaveRedlichKwong<Specie>& srk
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(srk);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        b_ = b_; 
    }
}


template<class Specie>
inline void Foam::soaveRedlichKwong<Specie>::operator*=(const scalar s)
{
     Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class Specie>
Foam::soaveRedlichKwong<Specie> Foam::operator+
(
    const soaveRedlichKwong<Specie>& srk1,
    const soaveRedlichKwong<Specie>& srk2
)
{
    Specie sp
    (
        static_cast<const Specie&>(srk1)
      + static_cast<const Specie&>(srk2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return soaveRedlichKwong<Specie>
        (
            sp,
            srk1.b_,
            srk1.coef1_,
            srk1.coef2_,
            srk1.coef3_
        );
    }
    else
    {
        const scalar b = srk1.b_;
        return soaveRedlichKwong<Specie>
        (
            sp,
            b,
            srk1.coef1_,
            srk1.coef2_,
            srk1.coef3_
        );
    }
}


template<class Specie>
Foam::soaveRedlichKwong<Specie> Foam::operator*
(
    const scalar s,
    const soaveRedlichKwong<Specie>& srk
)
{
    return soaveRedlichKwong<Specie>
    (
        s*static_cast<const Specie&>(srk),
        srk.b_,
        srk.coef1_,
        srk.coef2_,
        srk.coef3_
    );
}


template<class Specie>
Foam::soaveRedlichKwong<Specie> Foam::operator==
(
    const soaveRedlichKwong<Specie>& srk1,
    const soaveRedlichKwong<Specie>& srk2
)
{
    Specie sp
    (
        static_cast<const Specie&>(srk1)
     == static_cast<const Specie&>(srk2)
    );

    const scalar b = srk1.b_;
    return soaveRedlichKwong<Specie>
    (
        sp,
        b,
        srk1.coef1_,
        srk1.coef2_,
        srk1.coef3_
    );
}


// ************************************************************************* //
