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
#include "thermodynamicConstants.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::chungKineticTransport<Thermo>::calculate
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
inline Foam::scalar Foam::chungKineticTransport<Thermo>::mu
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
inline Foam::scalar Foam::chungKineticTransport<Thermo>::kappa
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
 inline Foam::scalar Foam::chungKineticTransport<Thermo>::alphah
 (
     const scalar p,
     const scalar T
 ) const
 {
 
     return kappa(p, T)/this->Cp(p, T);
 }



//- Return mass Diffusivity of specie ith [m^2/s] based on Fick's law using 
//  binary diffusion coefficients based on Standard Kinetic Theory
 template<class Thermo>
 inline Foam::scalar Foam::chungKineticTransport<Thermo>::Dimix
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
                       //(0.002628*sqrt(pow(T,3)/(2*Mij_[j][speciei])))/
                       (0.00266*sqrt(pow(T,3)/(2*Mij_[j][speciei])))/
                       ((p/1e5)*pow(sigmaij_[j][speciei],2)*
                        OMEGA11(T/epsilonijOverKb_[j][speciei], deltaij_[j][speciei]))                         
                   );
             }
         }  
         //Dimix = (1-Ymd_[speciei])/sum;
         if (sum == 0.0){Dimix = 0;} 
         else {Dimix = (1-Ymd_[speciei])/(sum);}
 
     return Dimix/1e4;
 }


// Quadratic interpolation
template<class Thermo>
inline Foam::scalar Foam::chungKineticTransport<Thermo>::qInterpolate
(
    scalar t,
    List<scalar> X,
    List<scalar> Y
) const
{
    scalar value = 0.0;
    scalar VAL1=0.0, VAL2=0.0, FAC1=0.0, FAC2=0.0;

    VAL1 = Y[0] + (t-X[0])*(Y[1]-Y[0]) / (X[1]-X[0]);
    VAL2 = Y[1] + (t-X[1])*(Y[2]-Y[1]) / (X[2]-X[1]);
    FAC1 = 0.5*(t-X[0]) / (X[1]-X[0]);
    FAC2 = 0.5*(X[2]-t) / (X[2]-X[1]);
   
    if ( t >= X[1])
    {
        value = (VAL1*FAC2+VAL2)/(1.0+FAC2);
    }
    else
    {
        value = (VAL2*FAC1+VAL1)/(1.0+FAC1);
    }
    return value; 
} 


//- Return collision integral from tablulated data
//- (Table IV, Monchick and Mason, Journal of Chem. Phys. 35: 1676 (1961))
template<class Thermo>
inline Foam::scalar Foam::chungKineticTransport<Thermo>::OMEGA11
(
    scalar TR,
    scalar DR
) const
{
    //table based on Stockmayer potentials
    // given in Table IV Monchick and Mason, J. Chem. Phys. 35:1676(1961)  
    const scalar Tstar[37] = 
    { 0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0, 
      1.2,  1.4,  1.6,  1.8,  2.0,  2.5,  3.0,  3.5,  4.0,  5.0,
      6.0,  7.0,  8.0,  9.0,  10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 
      25.0, 30.0, 35.0, 40.0, 50.0, 75.0, 100.0};  

    const scalar delta[8] = {0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5};

    const scalar table[37][8] = {    
       {4.008 , 4.002 , 4.655 , 5.52  , 6.454 , 8.214 , 9.824 , 11.31},
       {3.130 , 3.164 , 3.355 , 3.721 , 4.198 , 5.23  , 6.225 , 7.160},
       {2.649 , 2.657 , 2.77  , 3.002 , 3.319 , 4.054 , 4.785 , 5.483},
       {2.314 , 2.32  , 2.402 , 2.572 , 2.812 , 3.386 , 3.972 , 4.539},
       {2.066 , 2.073 , 2.14  , 2.278 , 2.472 , 2.946 , 3.437 , 3.918},
       {1.877 , 1.885 , 1.944 , 2.06  , 2.225 , 2.628 , 3.054 , 3.747},
       {1.729 , 1.738 , 1.79  , 1.893 , 2.036 , 2.388 , 2.763 , 3.137},
       {1.6122, 1.622 , 1.67  , 1.76  , 1.886 , 2.198 , 2.535 , 2.872},
       {1.517 , 1.527 , 1.572 , 1.653 , 1.765 , 2.044 , 2.35  , 2.657},
       {1.44  , 1.45  , 1.49  , 1.564 , 1.665 , 1.917 , 2.196 , 2.4780},

       {1.3204, 1.33  , 1.364 , 1.425 , 1.51  , 1.72  , 1.956 , 2.199},
       {1.234 , 1.24  , 1.272 , 1.324 , 1.394 , 1.573 , 1.777 , 1.99 },
       {1.168 , 1.176 , 1.202 , 1.246 , 1.306 , 1.46  , 1.64  , 1.827},
       {1.1166, 1.124 , 1.146 , 1.185 , 1.237 , 1.372 , 1.53  , 1.7  },
       {1.075 , 1.082 , 1.102 , 1.135 , 1.181 , 1.3   , 1.441 , 1.592},
       {1.0006, 1.005 , 1.02  , 1.046 , 1.08  , 1.17  , 1.278 , 1.397},
       {0.95  , 0.9538, 0.9656, 0.9852, 1.012 , 1.082 , 1.168 , 1.265},
       {0.9131, 0.9162, 0.9256, 0.9413, 0.9626, 1.019 , 1.09  , 1.17 },
       {0.8845, 0.8871, 0.8948, 0.9076, 0.9252, 0.972 , 1.03  , 1.098},
       {0.8428, 0.8446, 0.850 , 0.859 , 0.8716, 0.9053, 0.9483, 0.9984},

       {0.813 , 0.8142, 0.8183, 0.825 , 0.8344, 0.8598, 0.8927, 0.9316},
       {0.7898, 0.791 , 0.794 , 0.7993, 0.8066, 0.8265, 0.8526, 0.8836},
       {0.7711, 0.772 , 0.7745, 0.7788, 0.7846, 0.8007, 0.822 , 0.8474},
       {0.7555, 0.7562, 0.7584, 0.7619, 0.7667, 0.78  , 0.7976, 0.8189},
       {0.7422, 0.743 , 0.7446, 0.7475, 0.7515, 0.7627, 0.7776, 0.796 },
       {0.72022,0.7206, 0.722 , 0.7241, 0.7271, 0.7354, 0.7464, 0.76  },
       {0.7025, 0.703 , 0.704 , 0.7055, 0.7078, 0.7142, 0.7228, 0.7334},
       {0.68776, 0.688, 0.6888, 0.6901, 0.6919, 0.697 , 0.704 , 0.7125},
       {0.6751, 0.6753, 0.676 , 0.677 , 0.6785, 0.6827, 0.6884, 0.6955},
       {0.664 , 0.6642, 0.6648, 0.6657, 0.6669, 0.6704, 0.6752, 0.681 },

       {0.6414, 0.6415, 0.6418, 0.6425, 0.6433, 0.6457, 0.649 , 0.653 },
       {0.6235, 0.6236, 0.6239, 0.6243, 0.6249, 0.6267, 0.629 , 0.632 },
       {0.60882,0.6089, 0.6091, 0.6094, 0.61  , 0.6112, 0.613 , 0.6154},
       {0.5964, 0.5964, 0.5966, 0.597 , 0.5972, 0.5983, 0.600 , 0.6017},
       {0.5763, 0.5763, 0.5764, 0.5766, 0.5768, 0.5775, 0.5785, 0.58  },
       {0.5415, 0.5415, 0.5416, 0.5416, 0.5418, 0.542 , 0.5424, 0.543 },
       {0.518 , 0.518 , 0.5182, 0.5184, 0.5184, 0.5185, 0.5186, 0.5187} };
    // end of Table

    scalar omega = 0.0;
    if 
    ( 
        (DR < -1.0E-5) || (DR > 2.5) || (TR < 0.09) || (TR > 500) || 
        ( (DR > 1.0E-5) && (TR > 75) )
    )
    {
        Info << " <<<<<<<<<<< ERROR: Collision integral is undefined <<<<<<<<<<" << endl;
    } 
    else
    {
        if (TR > 75)
        {
            omega = 0.623 - 0.136E-2*TR + 0.346E-5*TR*TR -0.343E-8*TR*TR*TR;
        }
        else
        {
            int rowIndex = 0, columnIndex = 0;
            List<scalar> threePointsX(3); 
            List<scalar> threePointsY(3); 
            List<scalar> interpolatedRowData(3); 

            // Step 3: find the row index
            if (TR <= 0.2){rowIndex = 1;}
            else
            {
                for (int i =1; i < 37; i++)
                {
                    if ((TR > Tstar[i-1]) && (TR <= Tstar[i]))
                    { rowIndex = i;}
                }  
            }
            
            // Step 4
            if ( DR >= 1.0E-5)
            {
                // find the column index
                if (DR <= 0.25) {columnIndex = 1;}  
                else 
                { 
                    for (int j = 1; j < 7; j++)
                    {
                        if ((delta[j-1] < DR) && (DR <= delta[j]))
                        {columnIndex = j;} 
                    }
                }

                // do the interpolation
                // step 4.2.1. interpolate to build data for correct column of DR
                forAll(interpolatedRowData, i) 
                {
                    forAll(threePointsX, j)
                    {
                        threePointsX[j] = delta[columnIndex-1+j];
                        threePointsY[j] = table[rowIndex-1+i][columnIndex-1+j];
                    }
                    interpolatedRowData[i] = qInterpolate(DR, threePointsX, threePointsY);
                }
                // step 4.2.2. interpolate to find the data for correct row of TR
                forAll(threePointsX, i)
                {
                    threePointsX[i] = Tstar[rowIndex-1+i];
                }

                omega = qInterpolate(TR, threePointsX, interpolatedRowData);
            }
            else 
            {
                forAll(threePointsX, i)
                {
                    threePointsX[i] = Tstar[rowIndex-1+i];
                    interpolatedRowData[i] = table[rowIndex-1+i][0];      
                }
                omega = qInterpolate(TR, threePointsX, interpolatedRowData);
            } 
        }
    }
    return omega;
}



template<class Thermo>
inline void Foam::chungKineticTransport<Thermo>::updateTRANS
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
    const List<List<scalar>>& epsilonijOverKb,
    const List<List<scalar>>& deltaij,
    const List<List<scalar>>& Mij,
    const List<List<scalar>>& sigmaij
)
{
    sigmaM_          = sigmaM;
    epsilonkM_       = epsilonkM;
    MM_              = MM;
    VcM_             = VcM;
    TcM_             = TcM;
    omegaM_          = omegaM;
    miuiM_           = miuiM;
    kappaiM_         = kappaiM;
    Ymd_             = Ymd;
    Xmd_             = Xmd;
    epsilonijOverKb_ = epsilonijOverKb;
    deltaij_         = deltaij;
    Mij_             = Mij;
    sigmaij_         = sigmaij;
}


 template<class Thermo>
 inline Foam::scalar Foam::chungKineticTransport<Thermo>::sigmaM() const
 { 
     return sigmaM_;
 } 

 template<class Thermo>
 inline Foam::scalar Foam::chungKineticTransport<Thermo>::epsilonkM() const
 { 
     return epsilonkM_;
 } 

 template<class Thermo>
 inline Foam::scalar Foam::chungKineticTransport<Thermo>::MM() const
 { 
     return MM_;
 } 

 template<class Thermo>
 inline Foam::scalar Foam::chungKineticTransport<Thermo>::VcM() const
 { 
     return VcM_;
 } 

 template<class Thermo>
 inline Foam::scalar Foam::chungKineticTransport<Thermo>::TcM() const
 { 
     return TcM_;
 } 

 template<class Thermo>
 inline Foam::scalar Foam::chungKineticTransport<Thermo>::omegaM() const
 { 
     return omegaM_;
 } 

 template<class Thermo>
 inline Foam::scalar Foam::chungKineticTransport<Thermo>::miuiM() const
 { 
     return miuiM_;
 } 

 template<class Thermo>
 inline Foam::scalar Foam::chungKineticTransport<Thermo>::kappaiM() const
 { 
     return kappaiM_;
 } 



 // for diffusivity
 template<class Thermo>
 inline Foam::scalar Foam::chungKineticTransport<Thermo>::epsilonOverKb() const
 {
     return epsilonOverKb_;
 }

 template<class Thermo>
 inline Foam::scalar Foam::chungKineticTransport<Thermo>::sigma() const
 {
     return sigma_;
 }

 template<class Thermo>
 inline Foam::scalar Foam::chungKineticTransport<Thermo>::polar() const
 {
     return polar_;
 }

 template<class Thermo>
 inline Foam::scalar Foam::chungKineticTransport<Thermo>::Kb() const
 {
     return Kb_;
 }



 template<class Thermo>
 inline Foam::List<Foam::scalar> 
 Foam::chungKineticTransport<Thermo>::Ymd() const
 { 
     return Ymd_;
 } 

 template<class Thermo>
 inline Foam::List<Foam::scalar> 
 Foam::chungKineticTransport<Thermo>::Xmd() const
 { 
     return Xmd_;
 } 

 template<class Thermo>
 inline Foam::List<Foam::List<Foam::scalar>> 
 Foam::chungKineticTransport<Thermo>::epsilonijOverKb() const
 { 
     return epsilonijOverKb_;
 } 

 template<class Thermo>
 inline Foam::List<Foam::List<Foam::scalar>> 
 Foam::chungKineticTransport<Thermo>::deltaij() const
 { 
     return deltaij_;
 } 

 template<class Thermo>
 inline Foam::List<Foam::List<Foam::scalar>> 
 Foam::chungKineticTransport<Thermo>::Mij() const
 { 
     return Mij_;
 } 

 template<class Thermo>
 inline Foam::List<Foam::List<Foam::scalar>> 
 Foam::chungKineticTransport<Thermo>::sigmaij() const
 { 
     return sigmaij_;
 } 

// ************************************************************************* //
