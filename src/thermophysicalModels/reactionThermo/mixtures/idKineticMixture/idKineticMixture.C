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

#include "idKineticMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::idKineticMixture<ThermoType>::constructSpeciesData
(
    const dictionary& thermoDict
)
{
    forAll(species_, i)
    {
        speciesData_.set
        (
            i,
            new ThermoType(thermoDict.subDict(species_[i]))
        );
    }

    return speciesData_[0];
}


template<class ThermoType>
void Foam::idKineticMixture<ThermoType>::correctMassFractions()
{
    // Multiplication by 1.0 changes Yt patches to "calculated"
    volScalarField Yt("Yt", 1.0*Y_[0]);

    for (label n=1; n<Y_.size(); n++)
    {
        Yt += Y_[n];
    }

    if (mag(max(Yt).value()) < rootVSmall)
    {
        FatalErrorInFunction
            << "Sum of mass fractions is zero for species " << this->species()
            << exit(FatalError);
    }

    forAll(Y_, n)
    {
        Y_[n] /= Yt;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::idKineticMixture<ThermoType>::idKineticMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const HashPtrTable<ThermoType>& thermoData,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicSpecieMixture(thermoDict, specieNames, mesh, phaseName),
    speciesData_(species_.size()),
    mixture_("mixture", *thermoData[specieNames[0]]),
    mixtureVol_("volMixture", *thermoData[specieNames[0]]),
    //
    numberOfSpecies_(species_.size()),
    ListSpeciesName_(species_.size()),
     // for diffusivity
    ListW_(species_.size()),
    EPSILONijOVERKB_(species_.size()),
    DELTAij_(species_.size()),
    Mij_(species_.size()),
    SIGMAij_(species_.size())
    //
{
    forAll(species_, i)
    {
        speciesData_.set
        (
            i,
            new ThermoType(*thermoData[species_[i]])
        );
    }

    //
    forAll(ListSpeciesName_, i)
    {
        ListSpeciesName_[i] = speciesData_[i].name();
    }
    //

    correctMassFractions();
 
    // precalculation for kinetic theory model

    //for Kinetic model
    forAll(ListW_, i)
    { 
        ListW_[i]  = speciesData_[i].W();
    }

    List<scalar> nEPSILONijOVERKB(species_.size());
    List<scalar> nDELTAij(species_.size());
    List<scalar> nMij(species_.size());
    List<scalar> nSIGMAij(species_.size());

    // Dip = dipole moment 
    const scalar DipMin = 1e-20;

    forAll(EPSILONijOVERKB_, i)
    {
        forAll(nEPSILONijOVERKB, j)
        {
            nMij[j] = 
               1/(1/speciesData_[i].W() + 1/speciesData_[j].W());

            nEPSILONijOVERKB[j] = 
               sqrt(speciesData_[i].epsilonOverKb()*speciesData_[j].epsilonOverKb());

            nSIGMAij[j] = 
               0.5*(speciesData_[i].sigma() + speciesData_[j].sigma());

            // calculate coeficient xi
            scalar xi = 1.0; 
            if ((speciesData_[i].miui() < DipMin) && (speciesData_[j].miui() > DipMin)) 
            {
                // miui_j > DipMin > miui_i --> j is polar, i is nonpolar              
                nDELTAij[j] = 0;
                xi = 1.0 + 
                     speciesData_[i].polar()*pow(speciesData_[j].miui(), 2)* 
                     sqrt(speciesData_[j].epsilonOverKb()/speciesData_[i].epsilonOverKb())/
                     (4*pow(speciesData_[i].sigma(), 3)*
                            (
                             1e+19*speciesData_[j].Kb()*speciesData_[j].epsilonOverKb()*
                             pow(speciesData_[j].sigma(), 3)
                            )
                     );
            }
            else if ((speciesData_[i].miui() > DipMin) && (speciesData_[j].miui() < DipMin))
            {
                // miui_j < DipMin < miui_i --> i is polar, j is nonpolar
                nDELTAij[j] = 0;
                xi = 1.0 +
                     speciesData_[j].polar()*pow(speciesData_[i].miui(), 2)*
                     sqrt(speciesData_[i].epsilonOverKb()/speciesData_[j].epsilonOverKb())/
                     (4*pow(speciesData_[j].sigma(), 3)*   
                            (
                             1e+19*speciesData_[i].Kb()*speciesData_[i].epsilonOverKb()*
                             pow(speciesData_[i].sigma(), 3)
                            )
                     ); 
            }
            else 
            {            
                xi = 1.0;
                nDELTAij[j] =
                   0.5*(speciesData_[i].miui()*speciesData_[j].miui())/
                   (nEPSILONijOVERKB[j]*1e+19*speciesData_[j].Kb()*pow(nSIGMAij[j], 3));
            }
           
            nEPSILONijOVERKB[j] = nEPSILONijOVERKB[j]*pow(xi, 2);
            nSIGMAij[j] = nSIGMAij[j]*pow(xi, -1/6);
        }

        EPSILONijOVERKB_[i] = nEPSILONijOVERKB;
        Mij_[i]             = nMij;
        SIGMAij_[i]         = nSIGMAij;
        DELTAij_[i]         = nDELTAij;
    }

}


template<class ThermoType>
Foam::idKineticMixture<ThermoType>::idKineticMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicSpecieMixture
    (
        thermoDict,
        thermoDict.lookup("species"),
        mesh,
        phaseName
    ),
    speciesData_(species_.size()),
    mixture_("mixture", constructSpeciesData(thermoDict)),
    mixtureVol_("volMixture", speciesData_[0])
{
    correctMassFractions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::idKineticMixture<ThermoType>::cellMixture
(
    const label celli
) const
{
    mixture_ = Y_[0][celli]*speciesData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixture_ += Y_[n][celli]*speciesData_[n];
    }

   //- Update coefficients for diffusivity mixture
    //- List of secies mole and mass fraction 
    List<scalar> X(Y_.size()); 
    List<scalar> Y(Y_.size()); 
    scalar sumXb = 0.0;  
    forAll(X, i)
    {
        sumXb = sumXb + Y_[i][celli]/ListW_[i]; 
    }  
    if (sumXb == 0){ sumXb = 1e-30;} 

    forAll(X, i)
    {
        X[i] = (Y_[i][celli]/ListW_[i])/sumXb;
        Y[i] = Y_[i][celli];
        if(X[i] <= 0) { X[i] = 0; }
        if(Y[i] <= 0) { Y[i] = 0; }
    }

    scalar WmixCorrect = 0.0, sumXcorrected = 0.0;
    forAll(X, i)
    {
        X[i] = X[i] + 1e-40;
        sumXcorrected = sumXcorrected + X[i];
    }
    
    forAll(X, i)
    {
        X[i] = X[i]/sumXcorrected;
        WmixCorrect = WmixCorrect + X[i]*ListW_[i];
    }

    forAll(Y, i)
    {
        Y[i] = X[i]*ListW_[i]/WmixCorrect;
    }

    // Update coefficients for mixture of Kinetic model
    mixture_.updateTRANS(Y, X, EPSILONijOVERKB_, DELTAij_, Mij_, SIGMAij_);

    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::idKineticMixture<ThermoType>::patchFaceMixture
(
    const label patchi,
    const label facei
) const
{
    mixture_ = Y_[0].boundaryField()[patchi][facei]*speciesData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixture_ += Y_[n].boundaryField()[patchi][facei]*speciesData_[n];
    }

   //- Update coefficients for diffusivity mixture
    //- List of secies mole and mass fraction 
    List<scalar> X(Y_.size());
    List<scalar> Y(Y_.size());
    scalar sumXb = 0.0;
    forAll(X, i)
    {
        sumXb = sumXb + Y_[i].boundaryField()[patchi][facei]/ListW_[i];
    }
    if (sumXb == 0){ sumXb = 1e-30;}

    forAll(X, i)
    {
        X[i] = (Y_[i].boundaryField()[patchi][facei]/ListW_[i])/sumXb;
        Y[i] = Y_[i].boundaryField()[patchi][facei];
        if(X[i] <= 0) { X[i] = 0; } 
        if(Y[i] <= 0) { Y[i] = 0; }
    }

    scalar WmixCorrect = 0.0, sumXcorrected = 0.0;
    forAll(X, i)
    {
        X[i] = X[i] + 1e-40;
        sumXcorrected = sumXcorrected + X[i];
    }

    forAll(X, i)
    {
        X[i] = X[i]/sumXcorrected;
        WmixCorrect = WmixCorrect + X[i]*ListW_[i];
    }

    forAll(Y, i)
    {
        Y[i] = X[i]*ListW_[i]/WmixCorrect;
    }

    // Update coefficients for mixture of Kinetic model
    mixture_.updateTRANS(Y, X, EPSILONijOVERKB_, DELTAij_, Mij_, SIGMAij_);

    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::idKineticMixture<ThermoType>::cellVolMixture
(
    const scalar p,
    const scalar T,
    const label celli
) const
{
    scalar rhoInv = 0.0;
    forAll(speciesData_, i)
    {
        rhoInv += Y_[i][celli]/speciesData_[i].rho(p, T);
    }

    mixtureVol_ =
        Y_[0][celli]/speciesData_[0].rho(p, T)/rhoInv*speciesData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixtureVol_ +=
            Y_[n][celli]/speciesData_[n].rho(p, T)/rhoInv*speciesData_[n];
    }

    return mixtureVol_;
}


template<class ThermoType>
const ThermoType& Foam::idKineticMixture<ThermoType>::
patchFaceVolMixture
(
    const scalar p,
    const scalar T,
    const label patchi,
    const label facei
) const
{
    scalar rhoInv = 0.0;
    forAll(speciesData_, i)
    {
        rhoInv +=
            Y_[i].boundaryField()[patchi][facei]/speciesData_[i].rho(p, T);
    }

    mixtureVol_ =
        Y_[0].boundaryField()[patchi][facei]/speciesData_[0].rho(p, T)/rhoInv
      * speciesData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixtureVol_ +=
            Y_[n].boundaryField()[patchi][facei]/speciesData_[n].rho(p,T)
          / rhoInv*speciesData_[n];
    }

    return mixtureVol_;
}


template<class ThermoType>
void Foam::idKineticMixture<ThermoType>::read
(
    const dictionary& thermoDict
)
{
    forAll(species_, i)
    {
        speciesData_[i] = ThermoType(thermoDict.subDict(species_[i]));
    }
}


// ************************************************************************* //
