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
    *** 
\*---------------------------------------------------------------------------*/

#include "chungTakaTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
 
 //- Construct from components 
template<class Thermo>
inline Foam::chungTakaTransport<Thermo>::chungTakaTransport
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
)
:
    Thermo(t),
    sigmaM_(sigmaM),
    epsilonkM_(epsilonkM),
    MM_(MM),
    VcM_(VcM),
    TcM_(TcM),
    omegaM_(omegaM),
    miuiM_(miuiM),
    kappaiM_(kappaiM),
    // for diffusivity
    Ymd_(Ymd),
    Xmd_(Xmd),
    Tcmd_(Tcmd),
    Pcmd_(Pcmd),
    Mmd_(Mmd),
    sigmd_(sigmd)
{}


 //- Construct as named copy 
template<class Thermo>
inline Foam::chungTakaTransport<Thermo>::chungTakaTransport
(   
    const word& name,
    const chungTakaTransport& ch
)
:
    Thermo(name, ch),
    sigmaM_(ch.sigmaM_),
    epsilonkM_(ch.epsilonkM_),
    MM_(ch.MM_),
    VcM_(ch.VcM_),
    TcM_(ch.TcM_),
    omegaM_(ch.omegaM_),
    miuiM_(ch.miuiM_),
    kappaiM_(ch.kappaiM_),
    // for diffusivity
    Ymd_(ch.Ymd_),
    Xmd_(ch.Xmd_),
    Tcmd_(ch.Tcmd_),
    Pcmd_(ch.Pcmd_),
    Mmd_(ch.Mmd_),
    sigmd_(ch.sigmd_)
{}


 //- Construct and return a clone 
template<class Thermo>
inline Foam::autoPtr<Foam::chungTakaTransport<Thermo>>
Foam::chungTakaTransport<Thermo>::clone() const
{
    return autoPtr<chungTakaTransport<Thermo>>
    (
        new chungTakaTransport<Thermo>(*this)
    );
} 


 //- Selector from dictionary
template<class Thermo>
inline Foam::autoPtr<Foam::chungTakaTransport<Thermo>>
Foam::chungTakaTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<chungTakaTransport<Thermo>>
    (
        new chungTakaTransport<Thermo>(dict)
    );
}


 //- Construct from dictionary   
template<class Thermo>
Foam::chungTakaTransport<Thermo>::chungTakaTransport
(
    const dictionary& dict
)
:
    Thermo(dict),  
    sigmaM_(0.809*pow(this->Vc_, 1.0/3)),
    epsilonkM_(this->Tc_/1.2593), 
    MM_(this->molWeight_), 
    VcM_(this->Vc_), 
    TcM_(this->Tc_), 
    omegaM_(this->omega_), 
    miuiM_(this->miui_), 
    kappaiM_(this->kappai_), 
    // for diffusivity
    Ymd_(2),
    Xmd_(2),
    Tcmd_(2),
    Pcmd_(2),
    Mmd_(2),
    sigmd_(2)
{
  //- Temporary initialization
    // initialize Ymd_
    forAll(Ymd_, i) 
    {
        Ymd_[i] = this->Y();
    } 

    // initialize Xmd_
    forAll(Xmd_, i) 
    {
        Xmd_[i] = this->Y()/this->W();
    } 

    // initializes Tcmd_
    List<scalar> TcTemp(2);    
    forAll(TcTemp, i) 
    {
        TcTemp[i] = this->Tc_;
    } 

    forAll(Tcmd_, i) 
    {
        Tcmd_[i] = TcTemp;
    } 

    // initializes Pcmd_
    List<scalar> PcTemp(2);    
    forAll(PcTemp, i) 
    {
        PcTemp[i] = this->Pc_;
    } 

    forAll(Pcmd_, i) 
    {
        Pcmd_[i] = PcTemp;
    } 

    // initializes Mmd_
    List<scalar> MTemp(2);    
    forAll(MTemp, i) 
    {
        MTemp[i] = this->W();
    }

    forAll(Mmd_, i) 
    {
        Mmd_[i] = MTemp;
    } 

    // initializes sigmd_
    List<scalar> sigTemp(2);    
    forAll(sigTemp, i) 
    {
        sigTemp[i] = this->sigmvi_;
    } 

    forAll(sigmd_, i) 
    {
        sigmd_[i] = sigTemp;
    } 

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::chungTakaTransport<Thermo>::write(Ostream& os) const
{
    os  << this->rfSpecie::name() << endl
        << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.add("sigmaM", sigmaM_);
    dict.add("epsilonkM", epsilonkM_);
    dict.add("MM", MM_);
    dict.add("VcM", VcM_);
    dict.add("TcM", TcM_);
    dict.add("omegaM", omegaM_);
    dict.add("miuiM", miuiM_);
    dict.add("kappaiM", kappaiM_);
    
    os  << indent << dict.dictName() << dict
        << decrIndent << token::END_BLOCK << nl;
}



// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const chungTakaTransport<Thermo>& ch
)
{
    ch.write(os);
    return os;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::chungTakaTransport<Thermo>::operator=
(
    const chungTakaTransport<Thermo>& ch
)
{
    Thermo::operator=(ch);

    sigmaM_    = ch.sigmaM_;
    epsilonkM_ = ch.epsilonkM_;
    MM_        = ch.MM_;
    VcM_       = ch.VcM_;
    TcM_       = ch.TcM_;
    omegaM_    = ch.omegaM_;
    miuiM_     = ch.miuiM_;
    kappaiM_   = ch.kappaiM_;
    // for diffusivity
    Ymd_       = ch.Ymd_;
    Xmd_       = ch.Xmd_;
    Tcmd_      = ch.Tcmd_;
    Pcmd_      = ch.Pcmd_;
    Mmd_       = ch.Mmd_;
    sigmd_     = ch.sigmd_;
}


template<class Thermo>
inline void Foam::chungTakaTransport<Thermo>::operator+=
(
    const chungTakaTransport<Thermo>& ch
)
{
    Thermo::operator+=(ch);
    
        sigmaM_    = ch.sigmaM_;
        epsilonkM_ = ch.epsilonkM_;
        MM_        = ch.MM_;
        VcM_       = ch.VcM_;
        TcM_       = ch.TcM_;
        omegaM_    = ch.omegaM_;
        miuiM_     = ch.miuiM_;
        kappaiM_   = ch.kappaiM_;
        // for diffusivity 
        Ymd_       = ch.Ymd_; 
        Xmd_       = ch.Xmd_; 
        Tcmd_      = ch.Tcmd_; 
        Pcmd_      = ch.Pcmd_; 
        Mmd_       = ch.Mmd_; 
        sigmd_     = ch.sigmd_; 
}


template<class Thermo>
inline void Foam::chungTakaTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::chungTakaTransport<Thermo> Foam::operator+
(
    const chungTakaTransport<Thermo>& ch1,
    const chungTakaTransport<Thermo>& ch2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ch1) + static_cast<const Thermo&>(ch2)
    );

        return chungTakaTransport<Thermo>
        (
            t,
            ch1.sigmaM_,
            ch1.epsilonkM_,
            ch1.MM_,
            ch1.VcM_,
            ch1.TcM_,
            ch1.omegaM_,
            ch1.miuiM_,
            ch1.kappaiM_,
            // for diffusivity
            ch1.Ymd_,
            ch1.Xmd_,
            ch1.Tcmd_,
            ch1.Pcmd_,
            ch1.Mmd_,
            ch1.sigmd_
        );
}


template<class Thermo>
inline Foam::chungTakaTransport<Thermo> Foam::operator*
(
    const scalar s,
    const chungTakaTransport<Thermo>& ch
)
{
    return chungTakaTransport<Thermo>
    (
        s*static_cast<const Thermo&>(ch),
        ch.sigmaM_,
        ch.epsilonkM_,
        ch.MM_,
        ch.VcM_,
        ch.TcM_,
        ch.omegaM_,
        ch.miuiM_,
        ch.kappaiM_,
        // for diffusivity
        ch.Ymd_,
        ch.Xmd_,
        ch.Tcmd_,
        ch.Pcmd_,
        ch.Mmd_,
        ch.sigmd_
    );
}

// ************************************************************************* //
