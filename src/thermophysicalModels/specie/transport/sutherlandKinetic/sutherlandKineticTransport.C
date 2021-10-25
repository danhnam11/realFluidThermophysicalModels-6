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

#include "sutherlandKineticTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
Foam::scalar Foam::sutherlandKineticTransport<Thermo>::readCoeff
(
    const word& coeffName,
    const dictionary& dict
)
{
    return readScalar(dict.subDict("transport").lookup(coeffName));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- construct from components
template<class Thermo>
inline Foam::sutherlandKineticTransport<Thermo>::sutherlandKineticTransport
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
)
:
    Thermo(t),
    As_(As),
    Ts_(Ts),
    // for diffusivity
    epsilonOverKb_(epsilonOverKb),
    sigma_(sigma),
    miui_(miui),
    polar_(polar),
    Kb_(Kb),

    Ymd_(Ymd),
    Xmd_(Xmd),
    epsilonijOverKb_(epsilonijOverKb),
    deltaij_(deltaij),
    Mij_(Mij),
    sigmaij_(sigmaij)
{}

// construct as named copy
template<class Thermo>
inline Foam::sutherlandKineticTransport<Thermo>::sutherlandKineticTransport
(
    const word& name,
    const sutherlandKineticTransport& st
)
:
    Thermo(name, st),
    As_(st.As_),
    Ts_(st.Ts_),
    // for diffusivity
    epsilonOverKb_(st.epsilonOverKb_),
    sigma_(st.sigma_),
    miui_(st.miui_),
    polar_(st.polar_),
    Kb_(st.Kb_),    

    Ymd_(st.Ymd_),
    Xmd_(st.Xmd_),
    epsilonijOverKb_(st.epsilonijOverKb_),
    deltaij_(st.deltaij_),
    Mij_(st.Mij_),
    sigmaij_(st.sigmaij_)
{}


// construct from dictionary
template<class Thermo>
Foam::sutherlandKineticTransport<Thermo>::sutherlandKineticTransport
(
    const dictionary& dict
)
:
    Thermo(dict),
    As_(readCoeff("As", dict)),
    Ts_(readCoeff("Ts", dict)),
    // for diffusivity
    epsilonOverKb_(readScalar(dict.subDict("dataForKineticTrans").lookup("epsilonOverKb"))),
    sigma_(readScalar(dict.subDict("dataForKineticTrans").lookup("sigma"))),
    miui_(readScalar(dict.subDict("dataForKineticTrans").lookup("miui"))),
    polar_(readScalar(dict.subDict("dataForKineticTrans").lookup("alpha"))),
    Kb_(8.314510/(6.0221367*1E23)),    
 
    Ymd_(2),
    Xmd_(2),
    epsilonijOverKb_(2),
    deltaij_(2),
    Mij_(2),
    sigmaij_(2)
{

  //- Temporary initialization
    forAll(Ymd_, i) 
    {
        Ymd_[i] = this->Y();
        Xmd_[i] = this->Y()/this->W();
    } 

    forAll(epsilonijOverKb_, i) 
    {
        epsilonijOverKb_[i] = epsilonOverKb_;
        deltaij_[i] = 1.0;
        Mij_[i] = this->W();
        sigmaij_[i] = sigma_;
    } 

}

// construct and return a clone
template<class Thermo>
inline Foam::autoPtr<Foam::sutherlandKineticTransport<Thermo>>
Foam::sutherlandKineticTransport<Thermo>::clone() const
{
    return autoPtr<sutherlandKineticTransport<Thermo>>
    (
        new sutherlandKineticTransport<Thermo>(*this)
    );
}

// Selector from dictionary
template<class Thermo>
inline Foam::autoPtr<Foam::sutherlandKineticTransport<Thermo>>
Foam::sutherlandKineticTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<sutherlandKineticTransport<Thermo>>
    (
        new sutherlandKineticTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::sutherlandKineticTransport<Thermo>::write(Ostream& os) const
{
    os  << this->specie::name() << endl
        << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.add("As", As_);
    dict.add("Ts", Ts_);

    os  << indent << dict.dictName() << dict
        << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const sutherlandKineticTransport<Thermo>& st
)
{
    st.write(os);
    return os;
}



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::sutherlandKineticTransport<Thermo>::operator=
(
    const sutherlandKineticTransport<Thermo>& st
)
{
    Thermo::operator=(st);

    As_              = st.As_;
    Ts_              = st.Ts_;
    // for diffusivity
    epsilonOverKb_   = st.epsilonOverKb_;
    sigma_           = st.sigma_;
    miui_            = st.miui_;	
    polar_           = st.polar_;
    Kb_              = st.Kb_;

    Ymd_             = st.Ymd_;
    Xmd_             = st.Xmd_;
    epsilonijOverKb_ = st.epsilonijOverKb_;
    deltaij_         = st.deltaij_;
    Mij_             = st.Mij_;
    sigmaij_         = st.sigmaij_;
}


template<class Thermo>
inline void Foam::sutherlandKineticTransport<Thermo>::operator+=
(
    const sutherlandKineticTransport<Thermo>& st
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(st);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        scalar Y2 = st.Y()/this->Y();

        As_ = Y1*As_ + Y2*st.As_;
        Ts_ = Y1*Ts_ + Y2*st.Ts_;
    }
        // for diffusivity
        epsilonOverKb_   = st.epsilonOverKb_;
        sigma_           = st.sigma_;
        miui_            = st.miui_;
        polar_           = st.polar_;
        Kb_              = st.Kb_;

        Ymd_             = st.Ymd_; 
        Xmd_             = st.Xmd_; 
        epsilonijOverKb_ = st.epsilonijOverKb_; 
        deltaij_         = st.deltaij_; 
        Mij_             = st.Mij_; 
        sigmaij_         = st.sigmaij_;
}


template<class Thermo>
inline void Foam::sutherlandKineticTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}



// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::sutherlandKineticTransport<Thermo> Foam::operator+
(
    const sutherlandKineticTransport<Thermo>& st1,
    const sutherlandKineticTransport<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );

    if (mag(t.Y()) < small)
    {
        return sutherlandKineticTransport<Thermo>
        (
            t,
            st1.As_,
            st1.Ts_,
            // for diffusivity
            st1.epsilonOverKb_,
            st1.sigma_,
            st1.miui_,
            st1.polar_,
            st1.Kb_,

            st1.Ymd_,
            st1.Xmd_,
            st1.epsilonijOverKb_,
            st1.deltaij_,
            st1.Mij_,
            st1.sigmaij_
        );
    }
    else
    {
        scalar Y1 = st1.Y()/t.Y();
        scalar Y2 = st2.Y()/t.Y();

        return sutherlandKineticTransport<Thermo>
        (
            t,
            Y1*st1.As_ + Y2*st2.As_,
            Y1*st1.Ts_ + Y2*st2.Ts_,
            // for diffusivity
            st1.epsilonOverKb_,
            st1.sigma_,
            st1.miui_,
            st1.polar_,
            st1.Kb_,

            st1.Ymd_,
            st1.Xmd_,
            st1.epsilonijOverKb_,
            st1.deltaij_,
            st1.Mij_,
            st1.sigmaij_
        );
    }
}


template<class Thermo>
inline Foam::sutherlandKineticTransport<Thermo> Foam::operator*
(
    const scalar s,
    const sutherlandKineticTransport<Thermo>& st
)
{
    return sutherlandKineticTransport<Thermo>
    (
        s*static_cast<const Thermo&>(st),
        st.As_,
        st.Ts_,
        // for diffusivity
        st.epsilonOverKb_,
        st.sigma_,
        st.miui_,
        st.polar_,
        st.Kb_,

        st.Ymd_,
        st.Xmd_,
        st.epsilonijOverKb_,
        st.deltaij_,
        st.Mij_,
        st.sigmaij_
    );
}


// ************************************************************************* //
