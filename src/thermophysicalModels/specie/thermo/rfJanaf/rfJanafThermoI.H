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

#include "rfJanafThermo.H"
#include "specie.H"
#include "rfSpecie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::rfJanafThermo<EquationOfState>::rfJanafThermo
(
    const EquationOfState& st,
    const scalar Tlow,
    const scalar Thigh,
    const scalar Tcommon,
    const typename rfJanafThermo<EquationOfState>::coeffArray& highCpCoeffs,
    const typename rfJanafThermo<EquationOfState>::coeffArray& lowCpCoeffs,
    const bool convertCoeffs
)
:
    EquationOfState(st),
    Tlow_(Tlow),
    Thigh_(Thigh),
    Tcommon_(Tcommon)
{
    if (convertCoeffs)
    {
        for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
        {
            highCpCoeffs_[coefLabel] = highCpCoeffs[coefLabel]*this->R();
            lowCpCoeffs_[coefLabel] = lowCpCoeffs[coefLabel]*this->R();
        }
    }
    else
    {
        for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
        {
            highCpCoeffs_[coefLabel] = highCpCoeffs[coefLabel];
            lowCpCoeffs_[coefLabel] = lowCpCoeffs[coefLabel];
        }
    }
}


template<class EquationOfState>
inline const typename Foam::rfJanafThermo<EquationOfState>::coeffArray&
Foam::rfJanafThermo<EquationOfState>::coeffs
(
    const scalar T
) const
{
    if (T < Tcommon_)
    {
        return lowCpCoeffs_;
    }
    else
    {
        return highCpCoeffs_;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::rfJanafThermo<EquationOfState>::rfJanafThermo
(
    const word& name,
    const rfJanafThermo& jt
)
:
    EquationOfState(name, jt),
    Tlow_(jt.Tlow_),
    Thigh_(jt.Thigh_),
    Tcommon_(jt.Tcommon_)
{
    for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        highCpCoeffs_[coefLabel] = jt.highCpCoeffs_[coefLabel];
        lowCpCoeffs_[coefLabel] = jt.lowCpCoeffs_[coefLabel];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::rfJanafThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    if (T < Tlow_ || T > Thigh_)
    {
        WarningInFunction
            << "attempt to use rfJanafThermo<EquationOfState>"
               " out of temperature range "
            << Tlow_ << " -> " << Thigh_ << ";  T = " << T
            << endl;

        return min(max(T, Tlow_), Thigh_);
    }
    else
    {
        return T;
    }
}


template<class EquationOfState>
inline Foam::scalar Foam::rfJanafThermo<EquationOfState>::Tlow() const
{
    return Tlow_;
}


template<class EquationOfState>
inline Foam::scalar Foam::rfJanafThermo<EquationOfState>::Thigh() const
{
    return Thigh_;
}


template<class EquationOfState>
inline Foam::scalar Foam::rfJanafThermo<EquationOfState>::Tcommon() const
{
    return Tcommon_;
}


template<class EquationOfState>
inline const typename Foam::rfJanafThermo<EquationOfState>::coeffArray&
Foam::rfJanafThermo<EquationOfState>::highCpCoeffs() const
{
    return highCpCoeffs_;
}


template<class EquationOfState>
inline const typename Foam::rfJanafThermo<EquationOfState>::coeffArray&
Foam::rfJanafThermo<EquationOfState>::lowCpCoeffs() const
{
    return lowCpCoeffs_;
}


template<class EquationOfState>
inline Foam::scalar Foam::rfJanafThermo<EquationOfState>::Cp
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return
        ((((a[4]*T + a[3])*T + a[2])*T + a[1])*T + a[0])
      + EquationOfState::Cp(p, T);
}


//- Nam
template<class EquationOfState>
inline Foam::scalar Foam::rfJanafThermo<EquationOfState>::CpIdeal
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return
        ((((a[4]*T + a[3])*T + a[2])*T + a[1])*T + a[0]);
     
}
//-------


template<class EquationOfState>
inline Foam::scalar Foam::rfJanafThermo<EquationOfState>::Ha
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
/* 
	//JW
	Info << "Absolute Enthalpy Called. The value is " << 
	((((a[4]/5.0*T + a[3]/4.0)*T + a[2]/3.0)*T + a[1]/2.0)*T + a[0])*T + a[5] + EquationOfState::H(p, T) 
	<< endl;	
*/
	return
    (
        ((((a[4]/5.0*T + a[3]/4.0)*T + a[2]/3.0)*T + a[1]/2.0)*T + a[0])*T
      + a[5]
    ) + EquationOfState::H(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::rfJanafThermo<EquationOfState>::Hs
(
    const scalar p,
    const scalar T
) const
{
/*
	//JW
    Info << "Sensible Enthalpy Called. The value is " << Ha(p, T) - Hc() << endl; 
*/
    return Ha(p, T) - Hc();
}


template<class EquationOfState>
inline Foam::scalar Foam::rfJanafThermo<EquationOfState>::Hc() const
{
    const coeffArray& a = lowCpCoeffs_;
/*   
	//JW
	Info << "Chemical Enthalpy Called. The value is " <<  
    ((((a[4]/5.0*Tstd + a[3]/4.0)*Tstd + a[2]/3.0)*Tstd + a[1]/2.0)*Tstd + a[0])*Tstd + a[5]  
    << endl; 
*/
	return
    (
        (
            (((a[4]/5.0*Tstd + a[3]/4.0)*Tstd + a[2]/3.0)*Tstd + a[1]/2.0)*Tstd
          + a[0]
        )*Tstd + a[5]
    );
}


template<class EquationOfState>
inline Foam::scalar Foam::rfJanafThermo<EquationOfState>::S
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return
    (
        (((a[4]/4.0*T + a[3]/3.0)*T + a[2]/2.0)*T + a[1])*T + a[0]*log(T)
      + a[6]
    ) + EquationOfState::S(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::rfJanafThermo<EquationOfState>::dGdT
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return -((a[0] + a[5]/T)/T + a[1]/2 + T*(a[2]/3 + T*(a[3]/4 + T*a[4]/5)));
}


template<class EquationOfState>
inline Foam::scalar Foam::rfJanafThermo<EquationOfState>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return
        (((4*a[4]*T + 3*a[3])*T + 2*a[2])*T + a[1]);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::rfJanafThermo<EquationOfState>::operator+=
(
    const rfJanafThermo<EquationOfState>& jt
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(jt);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = jt.Y()/this->Y();

        Tlow_ = max(Tlow_, jt.Tlow_);
        Thigh_ = min(Thigh_, jt.Thigh_);

        if
        (
            rfJanafThermo<EquationOfState>::debug
         && notEqual(Tcommon_, jt.Tcommon_)
        )
        {
            FatalErrorInFunction
                << "Tcommon " << Tcommon_ << " for "
                << (this->name().size() ? this->name() : "others")
                << " != " << jt.Tcommon_ << " for "
                << (jt.name().size() ? jt.name() : "others")
                << exit(FatalError);
        }

        for
        (
            label coefLabel=0;
            coefLabel<rfJanafThermo<EquationOfState>::nCoeffs_;
            coefLabel++
        )
        {
            highCpCoeffs_[coefLabel] =
                Y1*highCpCoeffs_[coefLabel]
              + Y2*jt.highCpCoeffs_[coefLabel];

            lowCpCoeffs_[coefLabel] =
                Y1*lowCpCoeffs_[coefLabel]
              + Y2*jt.lowCpCoeffs_[coefLabel];
        }
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::rfJanafThermo<EquationOfState> Foam::operator+
(
    const rfJanafThermo<EquationOfState>& jt1,
    const rfJanafThermo<EquationOfState>& jt2
)
{
    EquationOfState eofs = jt1;
    eofs += jt2;

    if (mag(eofs.Y()) < small)
    {
        return rfJanafThermo<EquationOfState>
        (
            eofs,
            jt1.Tlow_,
            jt1.Thigh_,
            jt1.Tcommon_,
            jt1.highCpCoeffs_,
            jt1.lowCpCoeffs_
        );
    }
    else
    {
        const scalar Y1 = jt1.Y()/eofs.Y();
        const scalar Y2 = jt2.Y()/eofs.Y();

        typename rfJanafThermo<EquationOfState>::coeffArray highCpCoeffs;
        typename rfJanafThermo<EquationOfState>::coeffArray lowCpCoeffs;

        for
        (
            label coefLabel=0;
            coefLabel<rfJanafThermo<EquationOfState>::nCoeffs_;
            coefLabel++
        )
        {
            highCpCoeffs[coefLabel] =
                Y1*jt1.highCpCoeffs_[coefLabel]
              + Y2*jt2.highCpCoeffs_[coefLabel];

            lowCpCoeffs[coefLabel] =
                Y1*jt1.lowCpCoeffs_[coefLabel]
              + Y2*jt2.lowCpCoeffs_[coefLabel];
        }

        if
        (
            rfJanafThermo<EquationOfState>::debug
         && notEqual(jt1.Tcommon_, jt2.Tcommon_)
        )
        {
            FatalErrorInFunction
                << "Tcommon " << jt1.Tcommon_ << " for "
                << (jt1.name().size() ? jt1.name() : "others")
                << " != " << jt2.Tcommon_ << " for "
                << (jt2.name().size() ? jt2.name() : "others")
                << exit(FatalError);
        }

        return rfJanafThermo<EquationOfState>
        (
            eofs,
            max(jt1.Tlow_, jt2.Tlow_),
            min(jt1.Thigh_, jt2.Thigh_),
            jt1.Tcommon_,
            highCpCoeffs,
            lowCpCoeffs
        );
    }
}


template<class EquationOfState>
inline Foam::rfJanafThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const rfJanafThermo<EquationOfState>& jt
)
{
    return rfJanafThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(jt),
        jt.Tlow_,
        jt.Thigh_,
        jt.Tcommon_,
        jt.highCpCoeffs_,
        jt.lowCpCoeffs_
    );
}


template<class EquationOfState>
inline Foam::rfJanafThermo<EquationOfState> Foam::operator==
(
    const rfJanafThermo<EquationOfState>& jt1,
    const rfJanafThermo<EquationOfState>& jt2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(jt1)
     == static_cast<const EquationOfState&>(jt2)
    );

    const scalar Y1 = jt2.Y()/eofs.Y();
    const scalar Y2 = jt1.Y()/eofs.Y();

    typename rfJanafThermo<EquationOfState>::coeffArray highCpCoeffs;
    typename rfJanafThermo<EquationOfState>::coeffArray lowCpCoeffs;

    for
    (
        label coefLabel=0;
        coefLabel<rfJanafThermo<EquationOfState>::nCoeffs_;
        coefLabel++
    )
    {
        highCpCoeffs[coefLabel] =
            Y1*jt2.highCpCoeffs_[coefLabel]
          - Y2*jt1.highCpCoeffs_[coefLabel];

        lowCpCoeffs[coefLabel] =
            Y1*jt2.lowCpCoeffs_[coefLabel]
          - Y2*jt1.lowCpCoeffs_[coefLabel];
    }

    if
    (
        rfJanafThermo<EquationOfState>::debug
     && notEqual(jt2.Tcommon_, jt1.Tcommon_)
    )
    {
        FatalErrorInFunction
            << "Tcommon " << jt2.Tcommon_ << " for "
            << (jt2.name().size() ? jt2.name() : "others")
            << " != " << jt1.Tcommon_ << " for "
            << (jt1.name().size() ? jt1.name() : "others")
            << exit(FatalError);
    }

    return rfJanafThermo<EquationOfState>
    (
        eofs,
        max(jt2.Tlow_, jt1.Tlow_),
        min(jt2.Thigh_, jt1.Thigh_),
        jt2.Tcommon_,
        highCpCoeffs,
        lowCpCoeffs
    );
}


// ************************************************************************* //
