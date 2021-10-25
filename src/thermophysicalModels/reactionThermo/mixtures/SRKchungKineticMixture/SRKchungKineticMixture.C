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

#include "SRKchungKineticMixture.H"
#include "thermodynamicConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::SRKchungKineticMixture<ThermoType>::constructSpeciesData
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
void Foam::SRKchungKineticMixture<ThermoType>::correctMassFractions()
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


template<class ThermoType>
void Foam::SRKchungKineticMixture<ThermoType>::calculateRealGas
(
    List<scalar> X, 
    scalar& bM, 
    scalar& coef1, 
    scalar& coef2,
    scalar& coef3,
    scalar& sigmaM,
    scalar& epsilonkM,
    scalar& MM,
    scalar& VcM,
    scalar& TcM,
    scalar& omegaM, 
    scalar& miuiM,
    scalar& kappaiM
) const
{
    forAll(BM_, i)
    {
        //- For SRK
        bM = bM + X[i]*BM_[i];
    }
    if (bM == 0) {bM = 1e-16;}

    scalar sigma3M = 0, epsilonkM0 = 0, omegaM0 = 0, MM0 = 0, miuiM0 = 0;

    forAll(COEF1_, i)
    {
        forAll(COEF1_, j)
        {
            //- For SRK
            coef1      = coef1 + (X[i]*X[j])*COEF1_[i][j];
            coef2      = coef2 + (X[i]*X[j])*COEF2_[i][j];
            coef3      = coef3 + (X[i]*X[j])*COEF3_[i][j];
            //- For visc. and cond. in Chung's model 
            sigma3M    = sigma3M + (X[i]*X[j])*SIGMA3M_[i][j];
            epsilonkM0 = epsilonkM0 + (X[i]*X[j])*EPSILONKM0_[i][j];
            omegaM0    = omegaM0 + (X[i]*X[j])*OMEGAM0_[i][j];
            MM0        = MM0 + (X[i]*X[j])*MM0_[i][j];
            miuiM0     = miuiM0 + (X[i]*X[j])*MIUIM0_[i][j];
            kappaiM    = kappaiM + (X[i]*X[j])*KAPPAIM_[i][j];
        }
    }
    //- For visc. and cond. in Chung's model
    if (sigma3M ==0){ sigma3M = 1e-30;}
    sigmaM = pow(sigma3M, 1.0/3);
    //if (sigma3M ==0){ sigma3M = 1e-16; sigmaM = 1e-16; }

    epsilonkM = epsilonkM0/sigma3M;
    if (epsilonkM ==0) { epsilonkM = 1e-30;}

    VcM = pow(sigmaM/0.809, 3);
    if (VcM ==0) { VcM = 1e-30; }

    TcM = 1.2593*epsilonkM;

    omegaM = omegaM0/sigma3M;

    MM = pow(MM0/(epsilonkM*pow(sigmaM, 2)), 2);
    if (MM ==0) { MM = 1e-30; }

    miuiM = pow(miuiM0*sigma3M*epsilonkM, 1.0/4);
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class ThermoType>
Foam::SRKchungKineticMixture<ThermoType>::SRKchungKineticMixture
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
    //- 
//    realFluidMode_(thermoDict.lookupOrDefault("realFluidMode", false)), 
    numberOfSpecies_(species_.size()),
    ListSpeciesName_(species_.size()),
    ListW_(species_.size()),
    //- For SRK
    BM_(species_.size()),
    COEF1_(species_.size()),
    COEF2_(species_.size()),
    COEF3_(species_.size()),
    //- For visc. and cond. in Chung's model
    SIGMA3M_(species_.size()),
    EPSILONKM0_(species_.size()),
    OMEGAM0_(species_.size()),
    MM0_(species_.size()),
    MIUIM0_(species_.size()),
    KAPPAIM_(species_.size()),
    //- For Diffusivity in standard kintetic theory model
    EPSILONijOVERKB_(species_.size()),
    DELTAij_(species_.size()),
    Mij_(species_.size()),
    SIGMAij_(species_.size())
    //- end
{ 

    forAll(species_, i)
    {
        speciesData_.set
        (
            i,
            new ThermoType(*thermoData[species_[i]])
        );
    }

    forAll(ListSpeciesName_, i)
    {
        ListSpeciesName_[i] = speciesData_[i].name();
    }

    correctMassFractions();

    //- Real gas precalculation
    forAll(BM_, i)
    { 
        ListW_[i]  = speciesData_[i].W();
        BM_[i]     = 0.08664*RR*speciesData_[i].Tc()/speciesData_[i].Pc();
    }

    List<scalar> nCOEF1(species_.size());
    List<scalar> nCOEF2(species_.size());
    List<scalar> nCOEF3(species_.size());

    List<scalar> nSIGMA3M(species_.size());
    List<scalar> nEPSILONKM0(species_.size());
    List<scalar> nOMEGAM0(species_.size());
    List<scalar> nMM0(species_.size());
    List<scalar> nMIUIM0(species_.size());
    List<scalar> nKAPPAIM(species_.size());

    forAll(COEF1_, i)
    {
        forAll(nCOEF1, j)
        {
            //- For SRK
            nCOEF1[j] = 
                sqrt((0.42747*pow(RR*speciesData_[i].Tc(), 2)/speciesData_[i].Pc())
                    *(0.42747*pow(RR*speciesData_[j].Tc(), 2)/speciesData_[j].Pc()))
               *(1+(0.48508+1.5517*speciesData_[i].omega()-0.15613*pow(speciesData_[i].omega(), 2)))
               *(1+(0.48508+1.5517*speciesData_[j].omega()-0.15613*pow(speciesData_[j].omega(), 2)));

            nCOEF2[j] = 
                sqrt((0.42747*pow(RR*speciesData_[i].Tc(), 2)/speciesData_[i].Pc())
                    *(0.42747*pow(RR*speciesData_[j].Tc(), 2)/speciesData_[j].Pc()))
               *(
                  (
                    (1.0+(0.48508+1.5517*speciesData_[j].omega()-0.15613*pow(speciesData_[j].omega(), 2)))
                   *(0.48508+1.5517*speciesData_[i].omega()-0.15613*pow(speciesData_[i].omega(), 2))
                   /sqrt(speciesData_[i].Tc())
                 )
                +(
                   (1.0+(0.48508+1.5517*speciesData_[i].omega()-0.15613*pow(speciesData_[i].omega(), 2)))
                  *(0.48508+1.5517*speciesData_[j].omega()-0.15613*pow(speciesData_[j].omega(), 2))
                  /sqrt(speciesData_[j].Tc())
                 )
               );

            nCOEF3[j] = 
                sqrt((0.42747*pow(RR*speciesData_[i].Tc(), 2)/speciesData_[i].Pc())
                    *(0.42747*pow(RR*speciesData_[j].Tc(), 2)/speciesData_[j].Pc())) 
               *(0.48508+1.5517*speciesData_[i].omega()-0.15613*pow(speciesData_[i].omega(), 2))
               *(0.48508+1.5517*speciesData_[j].omega()-0.15613*pow(speciesData_[j].omega(), 2))
               /sqrt(speciesData_[i].Tc()*speciesData_[j].Tc());

            //- For visc. and cond. in Chung's model
            nSIGMA3M[j] = 
                pow(
                     (0.809*pow(speciesData_[i].Vc(), 1.0/3.0))*(0.809*pow(speciesData_[j].Vc(), 1.0/3.0)), 
                     3.0/2
                   );        

            nEPSILONKM0[j] = 
                pow(0.809*pow(speciesData_[i].Vc(),1.0/3)*0.809*pow(speciesData_[j].Vc(),1.0/3), 3.0/2)
               *pow((speciesData_[i].Tc()/1.2593)*(speciesData_[j].Tc()/1.2593), 1.0/2);

            nOMEGAM0[j] = 
                pow(0.809*pow(speciesData_[i].Vc(),1.0/3)*0.809*pow(speciesData_[j].Vc(),1.0/3), 3.0/2)
               *0.5*(speciesData_[i].omega()+ speciesData_[j].omega());  

            nMM0[j] = 
                (0.809*pow(speciesData_[i].Vc(),1.0/3))*(0.809*pow(speciesData_[j].Vc(),1.0/3))
               *pow((speciesData_[i].Tc()/1.2593)*(speciesData_[j].Tc()/1.2593), 1.0/2)
               *pow(2*speciesData_[i].W()*speciesData_[j].W()/(speciesData_[i].W()+speciesData_[j].W()),1.0/2); 

            nMIUIM0[j] = 
                pow(speciesData_[i].miui(), 2)*pow(speciesData_[j].miui(), 2)
               /(pow(0.809*pow(speciesData_[i].Vc(),1.0/3)*0.809*pow(speciesData_[j].Vc(),1.0/3), 3.0/2)
                *pow((speciesData_[i].Tc()/1.2593)*(speciesData_[j].Tc()/1.2593), 0.5) 
                );

            nKAPPAIM[j] = 
                pow(speciesData_[i].kappai()*speciesData_[j].kappai(), 1.0/2);
        }

        COEF1_[i] = nCOEF1;
        COEF2_[i] = nCOEF2;
        COEF3_[i] = nCOEF3;

        SIGMA3M_[i]    = nSIGMA3M;
        EPSILONKM0_[i] = nEPSILONKM0;
        OMEGAM0_[i]    = nOMEGAM0;
        MM0_[i]        = nMM0;
        MIUIM0_[i]     = nMIUIM0;
        KAPPAIM_[i]    = nKAPPAIM;
       //- end of real gas precalculation
    }


    // for Kinetic model
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
                     //speciesData_[i].polar()*speciesData_[j].miui()*       //original
                     speciesData_[i].polar()*pow(speciesData_[j].miui(), 2)* //CHEMKIN
                     sqrt(speciesData_[j].epsilonOverKb()/speciesData_[i].epsilonOverKb())/
                     (4*pow(speciesData_[i].sigma(), 3)*
                        //sqrt(  // original
                            (    // CHEMKIN
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
                     //speciesData_[j].polar()*speciesData_[i].miui()*       //original
                     speciesData_[j].polar()*pow(speciesData_[i].miui(), 2)* //CHEMKIN
                     sqrt(speciesData_[i].epsilonOverKb()/speciesData_[j].epsilonOverKb())/
                     (4*pow(speciesData_[j].sigma(), 3)*   
                        //sqrt(    // original
                            (      // CHEMKIN
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
Foam::SRKchungKineticMixture<ThermoType>::SRKchungKineticMixture
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
const ThermoType& Foam::SRKchungKineticMixture<ThermoType>::cellMixture
(
    const label celli
) const
{
    mixture_ = Y_[0][celli]*speciesData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixture_ += Y_[n][celli]*speciesData_[n];
    }

   //- Update coefficients for real gas mixture
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

    //- Using calculateRealGas function
    scalar bM = 0, coef1 = 0, coef2 = 0, coef3 = 0;
    scalar sigmaM = 0, epsilonkM =0, VcM = 0, TcM = 0, omegaM = 0, MM = 0, miuiM = 0, kappaiM = 0;

    calculateRealGas
    (
        X, bM, coef1, coef2, coef3,
        sigmaM, epsilonkM, MM, VcM, TcM, omegaM, miuiM, kappaiM
    );

    // Update coefficients for mixture in SRK
    mixture_.updateEoS(bM, coef1, coef2, coef3); 

    //- For mass diffusivity 
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
   
    // Update coefficients for mixture in Chung + Kinetic model
    mixture_.updateTRANS(sigmaM, epsilonkM, MM, VcM, TcM, omegaM, miuiM, kappaiM,
                       Y, X, EPSILONijOVERKB_, DELTAij_, Mij_, SIGMAij_);

    return mixture_;
}



template<class ThermoType>
const ThermoType& Foam::SRKchungKineticMixture<ThermoType>::patchFaceMixture
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

    //- Update coefficients for real gas mixture
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

    //- Using calculateRealGas function
    scalar bM = 0, coef1 = 0, coef2 = 0, coef3 = 0;
    scalar sigmaM = 0, epsilonkM =0, VcM = 0, TcM = 0, omegaM = 0, MM = 0, miuiM = 0, kappaiM = 0;

    calculateRealGas
    (
        X, bM, coef1, coef2, coef3,
        sigmaM, epsilonkM, MM, VcM, TcM, omegaM, miuiM, kappaiM
    );

    // Update coefficients for mixture in SRK
    mixture_.updateEoS(bM, coef1, coef2, coef3); 

    //- For mass diffusivity 
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
 
    // Update coefficients for mixture in Chung + kinetic model
    mixture_.updateTRANS(sigmaM, epsilonkM, MM, VcM, TcM, omegaM, miuiM, kappaiM,
                       Y, X, EPSILONijOVERKB_, DELTAij_, Mij_, SIGMAij_);
    return mixture_;
}



template<class ThermoType>
const ThermoType& Foam::SRKchungKineticMixture<ThermoType>::cellVolMixture
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
const ThermoType& Foam::SRKchungKineticMixture<ThermoType>::
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
void Foam::SRKchungKineticMixture<ThermoType>::read
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
