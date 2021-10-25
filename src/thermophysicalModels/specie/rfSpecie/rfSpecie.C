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
#include "constants.H"

/* * * * * * * * * * * * * * * public constants  * * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(rfSpecie, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rfSpecie::rfSpecie(const dictionary& dict)
:
    name_(dict.dictName()),
    Y_(dict.subDict("specie").lookupOrDefault("massFraction", 1.0)),
    molWeight_(readScalar(dict.subDict("specie").lookup("molWeight"))),
    Tc_(readScalar(dict.subDict("dataForRealGas").lookup("Tc"))),
    Pc_(1e6*readScalar(dict.subDict("dataForRealGas").lookup("Pc"))),
    Vc_(readScalar(dict.subDict("dataForRealGas").lookup("Vc"))), 
    omega_(readScalar(dict.subDict("dataForRealGas").lookup("omega"))),
    kappai_(readScalar(dict.subDict("dataForRealGas").lookup("kappai"))),
    miui_(readScalar(dict.subDict("dataForRealGas").lookup("miui"))),
    sigmvi_(readScalar(dict.subDict("dataForRealGas").lookup("sigmvi")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rfSpecie::write(Ostream& os) const
{
    dictionary dict("specie");
    if (Y_ != 1)
    {
        dict.add("massFraction", Y_);
    }
    dict.add("molWeight", molWeight_);
    dict.add("Tc", Tc_); 
    dict.add("Pc", Pc_); 
    dict.add("Vc", Vc_); 
    dict.add("omega", omega_); 
    dict.add("kappai", kappai_); 
    dict.add("miui", miui_); 
    dict.add("sigmvi", sigmvi_); 
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const rfSpecie& st)
{
    st.write(os);
    os.check("Ostream& operator<<(Ostream& os, const rfSpecie& st)");
    return os;
}


// ************************************************************************* //
