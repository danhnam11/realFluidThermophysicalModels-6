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

#include "rfSpecie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// construct from component with name 
inline rfSpecie::rfSpecie
(
    const word& name,
    const scalar Y,
    const scalar molWeight,
    const scalar Tc,
    const scalar Pc,
    const scalar Vc,
    const scalar omega,
    const scalar kappai,
    const scalar miui,
    const scalar sigmvi
)
:
    name_(name),
    Y_(Y),
    molWeight_(molWeight),
    Tc_(Tc),
    Pc_(Pc),
    Vc_(Vc),
    omega_(omega),
    kappai_(kappai),
    miui_(miui),
    sigmvi_(sigmvi)
{}


// construct from component without name
inline rfSpecie::rfSpecie
(
    const scalar Y,
    const scalar molWeight,
    const scalar Tc,
    const scalar Pc,
    const scalar Vc,
    const scalar omega,
    const scalar kappai,
    const scalar miui,
    const scalar sigmvi
)
:
    Y_(Y),
    molWeight_(molWeight),
    Tc_(Tc),
    Pc_(Pc),
    Vc_(Vc),
    omega_(omega),
    kappai_(kappai),
    miui_(miui),
    sigmvi_(sigmvi)
{}


// construct as named copy
inline rfSpecie::rfSpecie(const word& name, const rfSpecie& st)
:
    name_(name),
    Y_(st.Y_),
    molWeight_(st.molWeight_),
    Tc_(st.Tc_),
    Pc_(st.Pc_),
    Vc_(st.Vc_),
    omega_(st.omega_),
    kappai_(st.kappai_),
    miui_(st.miui_),
    sigmvi_(st.sigmvi_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const word& rfSpecie::name() const
{
    return name_;
}


inline scalar rfSpecie::W() const
{
    return molWeight_;
}


inline scalar rfSpecie::Y() const
{
    return Y_;
}


inline scalar rfSpecie::R() const
{
    return RR/molWeight_;
}


// Functions to be used in real-fluid evaluation
inline scalar rfSpecie::Tc() const
{
    return Tc_;
}

inline scalar rfSpecie::Pc() const
{
    return Pc_;
}

inline scalar rfSpecie::Vc() const
{
    return Vc_;
}

inline scalar rfSpecie::omega() const
{
    return omega_;
}

inline scalar rfSpecie::kappai() const
{
    return kappai_;
}

inline scalar rfSpecie::miui() const
{
    return miui_;
}

inline scalar rfSpecie::sigmvi() const
{
    return sigmvi_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void rfSpecie::operator=(const rfSpecie& st)
{
    // name_ = st.name_;
    Y_ = st.Y_;
    molWeight_ = st.molWeight_;
    Tc_ = st.Tc_;
    Pc_ = st.Pc_;
    Vc_ = st.Vc_;
    omega_ = st.omega_;
    kappai_ = st.kappai_;
    miui_ = st.miui_;
    sigmvi_ = st.sigmvi_;
}


inline void rfSpecie::operator+=(const rfSpecie& st)
{
    const scalar sumY = Y_ + st.Y_;
    if (mag(sumY) > small)
    {
        molWeight_ = sumY/(Y_/molWeight_ + st.Y_/st.molWeight_);
        Tc_ = st.Tc_;
        Pc_ = st.Pc_;
        Vc_ = st.Vc_;  
        omega_ = st.omega_; 
        kappai_ = st.kappai_; 
        miui_ = st.miui_; 
        sigmvi_ = st.sigmvi_; 
    }

    Y_ = sumY;
}


inline void rfSpecie::operator*=(const scalar s)
{
    Y_ *= s;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline rfSpecie operator+(const rfSpecie& st1, const rfSpecie& st2)
{
    const scalar sumY = max(st1.Y_ + st2.Y_, small);

    if (mag(sumY) > small)
    {
        return rfSpecie
        (
            sumY,
            sumY/(st1.Y_/st1.molWeight_ + st2.Y_/st2.molWeight_),
            (st1.Tc_*st1.Y_ + st2.Tc_*st2.Y_)/sumY,
            (st1.Pc_*st1.Y_ + st2.Pc_*st2.Y_)/sumY,
            (st1.Vc_*st1.Y_ + st2.Vc_*st2.Y_)/sumY,
            (st1.omega_*st1.Y_ + st2.omega_*st2.Y_)/sumY,
            (st1.kappai_*st1.Y_ + st2.kappai_*st2.Y_)/sumY,
            (st1.miui_*st1.Y_ + st2.miui_*st2.Y_)/sumY,
            (st1.sigmvi_*st1.Y_ + st2.sigmvi_*st2.Y_)/sumY
        );
    }
    else
    {
        return st1;
    }
}


inline rfSpecie operator*(const scalar s, const rfSpecie& st)
{
    return rfSpecie
    (
        s*st.Y_,
        st.molWeight_,
        st.Tc_,
        st.Pc_,
        st.Vc_,
        st.omega_,
        st.kappai_,
        st.miui_,
        st.sigmvi_
    );
}


inline rfSpecie operator==(const rfSpecie& st1, const rfSpecie& st2)
{
    scalar diffY = st2.Y_ - st1.Y_;
    if (mag(diffY) < small)
    {
        diffY = small;
    }

    const scalar diffRW =
        st2.Y_/st2.molWeight_ - st1.Y_/st1.molWeight_;

    scalar molWeight = great;
    scalar Tc = great;
    scalar Pc = great;
    scalar Vc = great;
    scalar omega = great;
    scalar kappai = great;
    scalar miui = great;
    scalar sigmvi = great;

    if (mag(diffRW) > small)
    {
        molWeight = diffY/diffRW;
        Tc = (st1.Tc_ + st2.Tc_)/2;
        Pc = (st1.Pc_ + st2.Pc_)/2;
        Vc = (st1.Vc_ + st2.Vc_)/2;
        omega = (st1.omega_ + st2.omega_)/2;
        kappai = (st1.kappai_ + st2.kappai_)/2;
        miui = (st1.miui_ + st2.miui_)/2;
        sigmvi = (st1.sigmvi_ + st2.sigmvi_)/2;
    }

    return rfSpecie
    (
        diffY,
        molWeight,
        Tc,
        Pc,
        Vc,
        omega,
        kappai,
        miui,
        sigmvi
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
