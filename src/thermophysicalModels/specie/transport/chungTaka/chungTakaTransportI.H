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

//#include "specie.H"
#include "rfSpecie.H"
#include "thermodynamicConstants.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::chungTakaTransport<Thermo>::calculate
(
    scalar p,
    scalar T,
    scalar& Tstar,
    scalar& mur,
    scalar& eta0,
    scalar& Y,
    scalar& G1
) const
{
    const scalar A=1.16145 , B=0.14874, C=0.52487 , D=0.77320 , E=2.16178; 
    const scalar F=2.43787 , G= -6.435e-4, H=7.27371, S=18.0323, W=-0.76830;

    // Calculate sigma by eq (4) 
//     scalar sigma = sigmaM_;
    // Calculate epsilon/k (=epsonk) in eq (5) 
//     scalar epsonk = epsilonkM_; 

    // Calculate T* (=Tstar) in eq (3) 
     Tstar = T/epsilonkM_; 

    // Calculate mu_r (=mur) in eq (8) 
     mur = 131.3*miuiM_/pow(VcM_*TcM_, 1.0/2.0); 
   
    // Calculate F_c (=Fc) in eq (7)
     scalar Fc = 1.0 - 0.2756*omegaM_ + 0.059035*pow(mur, 4.0) + kappaiM_; 

    // Calculate the reduced collision integral Omega* (=Omegast) in eq (2)  
     scalar Omegast = A*pow(Tstar, -B) + C*(exp(-D*Tstar)) + E*exp(-F*Tstar)
                          + G*pow(Tstar, B)*sin(S*pow(Tstar, W) - H); 
 
    // Calculate eta_0 (=eta0) in eq (6) 
     eta0 = (4.0785e-5)*pow(MM_*T, 1.0/2.0)*Fc/(Omegast*pow(VcM_, 2.0/3.0)); 

    // Call rho from based class and change unit from [kg/m3] to [mol/l] 
     scalar rho = 1e-3*(this->rho(p, T))/this->W(); 

    // Calculate Y term in eq (10) 
     Y =  rho*VcM_/6.0;

    // Calculate G1 term in eq (10) 
     G1 = (1.0 -  0.5*Y)/pow(1.0 - Y, 3.0);
}


//- Return Dynamic Viscosity [kg/m.s] based on Chung's model (1988) 
template<class Thermo>
inline Foam::scalar Foam::chungTakaTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    scalar sumX  = 0, finalMu = 0;
    forAll(Xmd_, i)
    {
         sumX = sumX + Xmd_[i];
    }

    if (sumX <= 1e-16)
    {
         //Info << "Sum of mole fraction = 0, then mu = 0" << endl;
         finalMu = 0;
  }
  else
  {
     scalar Tstar=0, mur=0, eta0=0, Y=0, G1=0;
     calculate(p, T, Tstar, mur, eta0, Y, G1); 
    
    // Specify the coefficients in table II, a0(i) -> a3(i), i=1,10 
     const scalar  a0[10] = {6.32402, 0.0012102, 5.28346, 6.62263, 19.74540, 
	                     -1.89992, 24.27450, 0.79716, -0.23816, 0.068629 };
      const scalar a1[10] = {50.41190, -0.0011536, 254.209, 38.0957, 7.63034, 
	                     -12.53670, 3.44945, 1.11764, 0.067695,  0.34793 };
      const scalar a2[10] = {-51.68, -0.0062571, -168.481, -8.46414, -14.35440, 
	                      4.98529, -11.29130, 0.012348, -0.8163, 0.59256 }; 
      const scalar a3[10] = {1189.02, 0.037283, 3898.270, 31.4178, 31.52670,
                               -18.1507, 69.3464, -4.1161, 4.02528, -0.72663 };  

    // Calculate the coefficients Ai in eq (11)    
     scalar MA[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};  
 
    for (int i=0; i<10; i++ )
    {
       // Notice that MA[0] = A1, A[1] = A2, ..., A[9] = A10 in eq (11) 
       MA[i] = a0[i] + a1[i]*omegaM_ + a2[i]*pow(mur, 4.0) + a3[i]*kappaiM_; 
    } 
    
    // Calulate G2 term in eq (10) 
     scalar G2 = (MA[0]*(1.0-exp(-MA[3]*Y))/Y + MA[1]*G1*exp(MA[4]*Y)+MA[2]*G1)/
                      (MA[0]*MA[3] + MA[1] + MA[2]); 
    
    // Calculate eta_k (=etak) in eq (10)
     scalar etak = eta0*(1.0/G2 + MA[5]*Y);  

    // Calculate eta_p (=etap) in eq (10)
     scalar etap = ((40.785e-6/sqrt(Tstar))*sqrt(MM_*T)/pow(VcM_, 2.0/3.0)) 
                        *MA[6]*Y*Y*G2*exp(MA[7]+MA[8]/Tstar + MA[9]/sqr(Tstar)); 
     // Return to final viscosity - mu (=eta) and change unit from poise [P] to [kg/ms] in eq (10)   
     finalMu = 0.1*(etak+etap); 
  }
    return finalMu; 
}


//- Return Thermal Conductivity [W/mK] based on Chung et al. model (1988) 
template<class Thermo>
inline Foam::scalar Foam::chungTakaTransport<Thermo>::kappa
(
    const scalar p, const scalar T
) const
{
    scalar sumX  = 0, kappa = 0;
    forAll(Xmd_, i)
    {
        sumX = sumX + Xmd_[i];
    }

    if (sumX <= 0.001)
    {
        kappa = 0;
  }
  else
  {
     scalar Tstar=0, mur=0, eta0=0, Y=0, G1=0;
     calculate(p, T, Tstar, mur, eta0, Y, G1); 

    // Specify the coefficients in table IV, b0(i) -> b3(i), i=1,7 
    const scalar b0[7] = {2.41657, -0.50924, 6.61069, 14.54250, 0.79274, 
                          -5.8634, 81.171}; 
    const scalar b1[7] = {0.74824, -1.50936, 5.62073, -8.91387, 0.82019, 
                          12.8005, 114.1580}; 
    const scalar b2[7] = {-0.91858, -49.99120, 64.7599, -5.63794, -0.69369, 
                           9.58926, -60.841}; 
    const scalar b3[7] = {121.721, 69.9834, 27.0389, 74.3435, 6.31734, 
                          -65.5292, 466.7750}; 

    // Calculate the coefficients Bi in eq (13) 
     scalar MB[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};  
 
    for (int j=0; j<7; j++ )
    {
       // Notice that MB[0] = B1, B[1] = B2, ..., B[6] = B7 in eq (13) 
       MB[j] = b0[j] + b1[j]*omegaM_ + b2[j]*pow(mur, 4.0) + b3[j]*kappaiM_; 
    } 
    
    // Call heat capacity at constant volume and change unit from [J/kgK] to [cal/mol.K]     
     scalar CvIdeal = this->CpIdeal(p,T) - this->R();
     scalar Cv = CvIdeal*MM_*0.2388e-3; 

    // Call universal gas constat from the paper [cal/mol.K] 
     //scalar R = 1.987; 
     scalar R = RR/(4.186*1e3);

    // Calculate alpha (=alpha) in eq (9) 
     scalar alpha = (Cv/R) - 3.0/2.0; 
     
    // Calculate beta (=beta) in eq (9) 
     scalar beta = 0.7862-0.7109*omegaM_ + 1.1368*pow(omegaM_, 2.0); 

    // Calculate Tr (=Tr) in eq (9)
     scalar Tr = T/TcM_; 

    // Calculate Z (=Z) in eq (9)
     scalar Z = 2.0 + 10.5*pow(Tr, 2.0); 

    // Calculate Psi (=Psi) in eq (9) 
     scalar Psi = 1.0 + alpha*(0.215 + 0.28288*alpha - 1.061*beta + 0.26665*Z)/ 
                          (0.6366 + beta*Z + 1.061*alpha*beta); 
  
    // Calculate lamda0 (=lamda0) in eq (9) 
     scalar lamda0 = 7.452*(eta0/MM_)*Psi; 

    // Calculate H2 (=H2) in eq (12)      
     scalar H2 = (MB[0]*(1.0-exp(-MB[3]*Y))/Y + MB[1]*G1*exp(MB[4]*Y)+MB[2]*G1)/
                      (MB[0]*MB[3] + MB[1] + MB[2]); 
   
    // Calculate lamda_k(=lamdak) in eq (12)  
     scalar lamdak = lamda0*(1.0/H2 + MB[5]*Y); 

    // Calculate lamda_p(=lamdap) in eq (12)
     scalar lamdap = (3.039e-4*pow(TcM_/MM_, 1.0/2.0)/pow(VcM_, 2.0/3.0))*MB[6]*pow(Y, 2.0)*H2*pow(Tr, 1.0/2.0); 

    // Return the final thermal conductivity and change unit from [cal/cm.s.K] to [w/mK] in eq (12) 
     kappa = 418.6798*(lamdak+lamdap);  
  }
    return kappa;
}



//- Return thermal diffusivity of enthalpy [kg/ms]-(original from sutherlandTransport class) 
 template<class Thermo>
 inline Foam::scalar Foam::chungTakaTransport<Thermo>::alphah
 (
     const scalar p,
     const scalar T
 ) const
 {
 
     return kappa(p, T)/this->Cp(p, T);
 }



//- Return mass Diffusivity of specie ith [m^2/s] based on Fick's law using 
//  binary Diffusion coefficients based on Fuller and Takahashi models
 template<class Thermo>
 inline Foam::scalar Foam::chungTakaTransport<Thermo>::Dimix
 (
     label speciei,
     const scalar p,
     const scalar T
 ) const
 {
         scalar sum = 0.0; 
         scalar Dimix = 0.0; 
         forAll(Xmd_, j) 
         {
             if (speciei != j)
             {
             sum = sum + Xmd_[j]/
                         (
                          phi(p/(Pcmd_[j][speciei]), T/(Tcmd_[j][speciei]))
                          *0.001*pow(T, 1.75)*pow(Mmd_[j][speciei], 0.5)
                          /((p/101325)*pow(sigmd_[j][speciei], 2))
                         );
             }
         }  
         //Dimix = (1-Ymd_[speciei])/sum;
         if (sum == 0.0){Dimix = 0;} 
         else {Dimix = (1-Ymd_[speciei])/(sum);}
         
     return Dimix/1e4;
 }


//- Return coefficient phi based on Takahashi model
template<class Thermo>
inline Foam::scalar Foam::chungTakaTransport<Thermo>::phi
(
    scalar PRR,
    scalar TRR
) const
{ 
    const scalar A01 = 0.038042, B01 = 1.52267; 
    const scalar A02 = 0.067433, B02 = 2.16794;
    const scalar A03 = 0.098317, B03 = 2.42910;
    const scalar A04 = 0.137610, B04 = 2.77605;
    const scalar A05 = 0.175081, B05 = 2.98256;
    const scalar A06 = 0.216376, B06 = 3.11384;
    const scalar A08 = 0.314051, B08 = 3.50264;
//    const scalar C08 = 0.000000, E08 = 0.00000; //unused
    const scalar A10 = 0.385736, B10 = 3.07773;
    const scalar C10 = 0.141211, E10 = 13.45454;
    const scalar A12 = 0.514553, B12 = 3.54744;
    const scalar C12 = 0.278407, E12 = 14.00000;
    const scalar A14 = 0.599184, B14 = 3.61216;
    const scalar C14 = 0.372683, E14 = 10.00900;
    const scalar A16 = 0.557725, B16 = 3.41882;
    const scalar C16 = 0.504894, E16 = 8.57519;
    const scalar A18 = 0.593007, B18 = 3.18415;
    const scalar C18 = 0.678469, E18 = 10.37483;
    const scalar A20 = 0.696001, B20 = 3.37660;
    const scalar C20 = 0.665702, E20 = 11.21674;
    const scalar A25 = 0.790770, B25 = 3.27984;
    const scalar C25 = 0.000000, E25 = 0.00000;
    const scalar A30 = 0.502100, B30 = 2.39031;
    const scalar C30 = 0.602907, E30 = 6.19043;
    const scalar A40 = 0.837452, B40 = 3.23513;
    const scalar C40 = 0.000000, E40 = 0.00000;
    const scalar A50 = 0.890390, B50 = 3.13001;

    scalar Y=0, YA=0, YB=0; 
      TRR = (TRR < 1.0) ? 1.0 : TRR; 
      TRR = (TRR > 5.0) ? 5.0 : TRR; 
      PRR = (PRR > 5.0) ? 5.0 : PRR; 

      if(PRR < 0.1)
          { 
            Y = 1.01 * (1.0 - A01 * pow(TRR, -B01));
            Y = 1.0 + (Y - 1.0)/(0.1-0.0)*(PRR-0.0);
          } 
      else if((PRR >= 0.1) && (PRR < 0.2))
	  { 
            YA = 1.01 * (1.0 - A01 * pow(TRR, -B01)); 
            YB = 1.01 * (1.0 - A02 * pow(TRR, -B02)); 
            Y  = YA + (YB-YA)/(0.2-0.1)*(PRR-0.1); 
	  }
      else if((PRR >= 0.2) && (PRR < 0.3))
	  { 
            YA = 1.01 * (1.0 - A02 * pow(TRR, -B02)); 
            YB = 1.01 * (1.0 - A03 * pow(TRR, -B03)); 
            Y  = YA + (YB-YA)/(0.3-0.2)*(PRR-0.2); 
	  } 
      else if((PRR >= 0.3) && (PRR < 0.4))
	  {
            YA = 1.01 * (1.0 - A03 * pow(TRR, -B03)); 
            YB = 1.01 * (1.0 - A04 * pow(TRR, -B04));
            Y  = YA + (YB-YA)/(0.4-0.3)*(PRR-0.3); 
	  }
      else if((PRR >= 0.4) && (PRR < 0.5))
	  {
            YA = 1.01 * (1.0 - A04 * pow(TRR, -B04));
            YB = 1.01 * (1.0 - A05 * pow(TRR, -B05));
            Y  = YA + (YB-YA)/(0.5-0.4)*(PRR-0.4);
	  }
      else if((PRR >= 0.5) && (PRR < 0.6))
	  {
            YA = 1.01 * (1.0 - A05 * pow(TRR, -B05)); 
            YB = 1.01 * (1.0 - A06 * pow(TRR, -B06)); 
            Y  = YA + (YB-YA)/(0.6-0.5)*(PRR-0.5); 
	  }
      else if((PRR >= 0.6) && (PRR < 0.8))
	  {
            YA = 1.01 * (1.0 - A06 * pow(TRR, -B06)); 
            YB = 1.01 * (1.0 - A08 * pow(TRR, -B08)); 
            Y  = YA + (YB-YA)/(0.8-0.6)*(PRR-0.6); 
          }
      else if((PRR >= 0.8) && (PRR < 1.0))
	  {
            YA = 1.01 * (1.0 - A08 * pow(TRR, -B08)); 
            YB = 1.02 * (1.0 - A10 * pow(TRR, -B10))*(1.0 - C10*pow(TRR, -E10)); 
            Y  = YA + (YB-YA)/(1.0-0.8)*(PRR-0.8);
	  }
      else if((PRR >= 1.0) && (PRR < 1.2))
	  { 
            YA = 1.02 * (1.0 - A10 * pow(TRR, -B10))*(1.0 - C10*pow(TRR, -E10));
            YB = 1.02 * (1.0 - A12 * pow(TRR, -B12))*(1.0 - C12*pow(TRR, -E12));
            Y  = YA + (YB-YA)/(1.2-1.0)*(PRR-1.0); 
	  }
      else if((PRR >= 1.2) && (PRR < 1.4))
	  {
            YA = 1.02 * (1.0 - A12 * pow(TRR, -B12))*(1.0 - C12*pow(TRR, -E12));
            YB = 1.02 * (1.0 - A14 * pow(TRR, -B14))*(1.0 - C14*pow(TRR, -E14));
            Y  = YA + (YB-YA)/(1.4-1.2)*(PRR-1.2); 
	  }
      else if((PRR >= 1.4) && (PRR < 1.6))
	  {
            YA = 1.02 * (1.0 - A14 * pow(TRR, -B14))*(1.0 - C14*pow(TRR, -E14)); 
            YB = 1.02 * (1.0 - A16 * pow(TRR, -B16))*(1.0 - C16*pow(TRR, -E16));
            Y  = YA + (YB-YA)/(1.6-1.4)*(PRR-1.4);
	  }
      else if((PRR >= 1.6) && (PRR < 1.8))
	  {
            YA = 1.02 * (1.0 - A16 * pow(TRR, -B16))*(1.0 - C16*pow(TRR, -E16));
            YB = 1.03 * (1.0 - A18 * pow(TRR, -B18))*(1.0 - C18*pow(TRR, -E18));
            Y  = YA + (YB-YA)/(1.8-1.6)*(PRR-1.6);
	  }
      else if((PRR >= 1.8) && (PRR < 2.0))
	  {
            YA = 1.03 * (1.0 - A18 * pow(TRR, -B18))*(1.0 - C18*pow(TRR, -E18));
            YB = 1.03 * (1.0 - A20 * pow(TRR, -B20))*(1.0 - C20*pow(TRR, -E20));
            Y  = YA + (YB-YA)/(2.0-1.8)*(PRR-1.8);
	  }
      else if((PRR >= 2.0) && (PRR < 2.5))
	  {
            YA = 1.03 * (1.0 - A20 * pow(TRR, -B20))*(1.0 - C20*pow(TRR, -E20));
            YB = 1.04 * (1.0 - A25 * pow(TRR, -B25))*(1.0 - C25*pow(TRR, -E25));
            Y  = YA + (YB-YA)/(2.5-2.0)*(PRR-2.0);
	  }
      else if((PRR >= 2.5) && (PRR < 3.0))
	  {
            YA = 1.04 * (1.0 - A25 * pow(TRR, -B25))*(1.0 - C25*pow(TRR, -E25));
            YB = 1.05 * (1.0 - A30 * pow(TRR, -B30))*(1.0 - C30*pow(TRR, -E30));
            Y  = YA + (YB-YA)/(3.0-2.5)*(PRR-2.5);
	  }
      else if((PRR >= 3.0) && (PRR < 4.0))
	  {
            YA = 1.05 * (1.0 - A30 * pow(TRR, -B30))*(1.0 - C30*pow(TRR, -E30));
            YB = 1.06 * (1.0 - A40 * pow(TRR, -B40))*(1.0 - C40*pow(TRR, -E40));
            Y  = YA + (YB-YA)/(4.0-3.0)*(PRR-3.0);
	  }
      else if((PRR >= 4.0) && (PRR < 5.0))
	  {
            YA = 1.06 * (1.0 - A40 * pow(TRR, -B40))*(1.0 - C40*pow(TRR, -E40));
            YB = 1.07 * (1.0 - A50 * pow(TRR, -B50));
			Y  = YA + (YB-YA)/(5.0-4.0)*(PRR-4.0);
	  }
      else if(PRR >= 5.0)
	  {
             Y = 1.07 * (1.0 - A50 * pow(TRR, -B50));
	  }

      if (Y ==0) { Y = 1e-16;}
      return Y;
} 


template<class Thermo>
inline void Foam::chungTakaTransport<Thermo>::updateTRANS
(
    const scalar& sigmaM,
    const scalar& epsilonkM,
    const scalar& MM,
    const scalar& VcM,
    const scalar& TcM,
    const scalar& omegaM,
    const scalar& miuiM,
    const scalar& kappaiM,
    const List<scalar>& Ymd, 
    const List<scalar>& Xmd, 
    const List<List<scalar>>& Tcmd,
    const List<List<scalar>>& Pcmd,
    const List<List<scalar>>& Mmd,
    const List<List<scalar>>& sigmd
)
{
    sigmaM_    = sigmaM;
    epsilonkM_ = epsilonkM;
    MM_        = MM;
    VcM_       = VcM;
    TcM_       = TcM;
    omegaM_    = omegaM;
    miuiM_     = miuiM;
    kappaiM_   = kappaiM;
    Ymd_       = Ymd;
    Xmd_       = Xmd;
    Tcmd_      = Tcmd;
    Pcmd_      = Pcmd;
    Mmd_       = Mmd;
    sigmd_     = sigmd;
}


 template<class Thermo>
 inline Foam::scalar Foam::chungTakaTransport<Thermo>::sigmaM() const
 { 
     return sigmaM_;
 } 

 template<class Thermo>
 inline Foam::scalar Foam::chungTakaTransport<Thermo>::epsilonkM() const
 { 
     return epsilonkM_;
 } 

 template<class Thermo>
 inline Foam::scalar Foam::chungTakaTransport<Thermo>::MM() const
 { 
     return MM_;
 } 

 template<class Thermo>
 inline Foam::scalar Foam::chungTakaTransport<Thermo>::VcM() const
 { 
     return VcM_;
 } 

 template<class Thermo>
 inline Foam::scalar Foam::chungTakaTransport<Thermo>::TcM() const
 { 
     return TcM_;
 } 

 template<class Thermo>
 inline Foam::scalar Foam::chungTakaTransport<Thermo>::omegaM() const
 { 
     return omegaM_;
 } 

 template<class Thermo>
 inline Foam::scalar Foam::chungTakaTransport<Thermo>::miuiM() const
 { 
     return miuiM_;
 } 

 template<class Thermo>
 inline Foam::scalar Foam::chungTakaTransport<Thermo>::kappaiM() const
 { 
     return kappaiM_;
 } 

 // for diffusivity
 template<class Thermo>
 inline Foam::List<Foam::scalar> Foam::chungTakaTransport<Thermo>::Ymd() const
 { 
     return Ymd_;
 } 

 template<class Thermo>
 inline Foam::List<Foam::scalar> Foam::chungTakaTransport<Thermo>::Xmd() const
 { 
     return Xmd_;
 } 

 template<class Thermo>
 inline Foam::List<Foam::List<Foam::scalar>> Foam::chungTakaTransport<Thermo>::Tcmd() const
 { 
     return Tcmd_;
 } 

 template<class Thermo>
 inline Foam::List<Foam::List<Foam::scalar>> Foam::chungTakaTransport<Thermo>::Pcmd() const
 { 
     return Pcmd_;
 } 

 template<class Thermo>
 inline Foam::List<Foam::List<Foam::scalar>> Foam::chungTakaTransport<Thermo>::Mmd() const
 { 
     return Mmd_;
 } 

 template<class Thermo>
 inline Foam::List<Foam::List<Foam::scalar>> Foam::chungTakaTransport<Thermo>::sigmd() const
 { 
     return sigmd_;
 } 

// ************************************************************************* //
