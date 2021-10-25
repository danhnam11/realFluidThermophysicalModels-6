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

#include "chungKineticTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
 
 //- Construct from components 
template<class Thermo>
inline Foam::chungKineticTransport<Thermo>::chungKineticTransport
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
    epsilonOverKb_(epsilonOverKb),
    sigma_(sigma),
    polar_(polar),
    Kb_(Kb),

    Ymd_(Ymd),
    Xmd_(Xmd),
    epsilonijOverKb_(epsilonijOverKb),
    deltaij_(deltaij),
    Mij_(Mij),
    sigmaij_(sigmaij)
{}


 //- Construct as named copy 
template<class Thermo>
inline Foam::chungKineticTransport<Thermo>::chungKineticTransport
(   
    const word& name,
    const chungKineticTransport& ch
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
    epsilonOverKb_(ch.epsilonOverKb_),
    sigma_(ch.sigma_),
    polar_(ch.polar_),
    Kb_(ch.Kb_),    

    Ymd_(ch.Ymd_),
    Xmd_(ch.Xmd_),
    epsilonijOverKb_(ch.epsilonijOverKb_),
    deltaij_(ch.deltaij_),
    Mij_(ch.Mij_),
    sigmaij_(ch.sigmaij_)
{}


 //- Construct and return a clone 
template<class Thermo>
inline Foam::autoPtr<Foam::chungKineticTransport<Thermo>>
Foam::chungKineticTransport<Thermo>::clone() const
{
    return autoPtr<chungKineticTransport<Thermo>>
    (
        new chungKineticTransport<Thermo>(*this)
    );
} 


 //- Selector from dictionary
template<class Thermo>
inline Foam::autoPtr<Foam::chungKineticTransport<Thermo>>
Foam::chungKineticTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<chungKineticTransport<Thermo>>
    (
        new chungKineticTransport<Thermo>(dict)
    );
}


 //- Construct from dictionary   
template<class Thermo>
Foam::chungKineticTransport<Thermo>::chungKineticTransport
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
    epsilonOverKb_(readScalar(dict.subDict("dataForKineticTrans").lookup("epsilonOverKb"))),
    sigma_(readScalar(dict.subDict("dataForKineticTrans").lookup("sigma"))),
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
        epsilonijOverKb_[i] = this->epsilonOverKb_;
        deltaij_[i] = 1.0;
        Mij_[i] = this->W();
        sigmaij_[i] = this->sigmvi_;
    } 

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::chungKineticTransport<Thermo>::write(Ostream& os) const
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
    const chungKineticTransport<Thermo>& ch
)
{
    ch.write(os);
    return os;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::chungKineticTransport<Thermo>::operator=
(
    const chungKineticTransport<Thermo>& ch
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
    epsilonOverKb_   = ch.epsilonOverKb_;
    sigma_           = ch.sigma_;
    polar_           = ch.polar_;
    Kb_              = ch.Kb_;

    Ymd_             = ch.Ymd_;
    Xmd_             = ch.Xmd_;
    epsilonijOverKb_ = ch.epsilonijOverKb_;
    deltaij_         = ch.deltaij_;
    Mij_             = ch.Mij_;
    sigmaij_         = ch.sigmaij_;
}


template<class Thermo>
inline void Foam::chungKineticTransport<Thermo>::operator+=
(
    const chungKineticTransport<Thermo>& ch
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
        epsilonOverKb_   = ch.epsilonOverKb_;
        sigma_           = ch.sigma_;
        polar_           = ch.polar_;
        Kb_              = ch.Kb_;

        Ymd_             = ch.Ymd_; 
        Xmd_             = ch.Xmd_; 
        epsilonijOverKb_ = ch.epsilonijOverKb_; 
        deltaij_         = ch.deltaij_; 
        Mij_             = ch.Mij_; 
        sigmaij_         = ch.sigmaij_; 
}


template<class Thermo>
inline void Foam::chungKineticTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::chungKineticTransport<Thermo> Foam::operator+
(
    const chungKineticTransport<Thermo>& ch1,
    const chungKineticTransport<Thermo>& ch2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ch1) + static_cast<const Thermo&>(ch2)
    );

        return chungKineticTransport<Thermo>
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
            ch1.epsilonOverKb_,
            ch1.sigma_,
            ch1.polar_,
            ch1.Kb_,

            ch1.Ymd_,
            ch1.Xmd_,
            ch1.epsilonijOverKb_,
            ch1.deltaij_,
            ch1.Mij_,
            ch1.sigmaij_
        );
}


template<class Thermo>
inline Foam::chungKineticTransport<Thermo> Foam::operator*
(
    const scalar s,
    const chungKineticTransport<Thermo>& ch
)
{
    return chungKineticTransport<Thermo>
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
        ch.epsilonOverKb_,
        ch.sigma_,
        ch.polar_,
        ch.Kb_,

        ch.Ymd_,
        ch.Xmd_,
        ch.epsilonijOverKb_,
        ch.deltaij_,
        ch.Mij_,
        ch.sigmaij_
    );
}

// ************************************************************************* //
