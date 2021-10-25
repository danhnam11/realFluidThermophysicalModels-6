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

#include "specie.H"
#include "thermodynamicConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
inline void Foam::sutherlandKineticTransport<Thermo>::calcCoeffs
(
    const scalar mu1, const scalar T1,
    const scalar mu2, const scalar T2
)
{
    scalar rootT1 = sqrt(T1);
    scalar mu1rootT2 = mu1*sqrt(T2);
    scalar mu2rootT1 = mu2*rootT1;

    Ts_ = (mu2rootT1 - mu1rootT2)/(mu1rootT2/T1 - mu2rootT1/T2);

    As_ = mu1*(1.0 + Ts_/T1)/rootT1;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::sutherlandKineticTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    return As_*::sqrt(T)/(1.0 + Ts_/T);
}


template<class Thermo>
inline Foam::scalar Foam::sutherlandKineticTransport<Thermo>::kappa
(
    const scalar p, const scalar T
) const
{
    scalar Cv_ = this->Cv(p, T);
    return mu(p, T)*Cv_*(1.32 + 1.77*this->R()/Cv_);
}


template<class Thermo>
inline Foam::scalar Foam::sutherlandKineticTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{

    return kappa(p, T)/this->Cp(p, T);
}


// for Diffusivity 
template<class Thermo>
inline Foam::scalar Foam::sutherlandKineticTransport<Thermo>::Dimix
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
                       (0.00266*sqrt(pow(T,3)/(2*Mij_[j][speciei])))/
                       ((p/1e5)*pow(sigmaij_[j][speciei],2)*
                        OMEGA11(T/epsilonijOverKb_[j][speciei], deltaij_[j][speciei]))                         
                   );
             }
         }  
         if (sum == 0.0){Dimix = 0;} 
         else {Dimix = (1-Ymd_[speciei])/(sum);}

     return Dimix/1e4;
} 


// Quadratic interpolation
template<class Thermo>
inline Foam::scalar Foam::sutherlandKineticTransport<Thermo>::qInterpolate
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
inline Foam::scalar Foam::sutherlandKineticTransport<Thermo>::OMEGA11
(
    scalar TR,
    scalar DR
) const
{
     //table based on Stockmayer potentials 
     //given in Table IV Monchick and Mason, J. Chem. Phys. 35:1676(1961)  
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
     //end of Table

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


// update coefficients for Diffusivity 
template<class Thermo>
inline void Foam::sutherlandKineticTransport<Thermo>::updateTRANS
(
    const List<scalar>& Ymd,
    const List<scalar>& Xmd,
    const List<List<scalar>>& epsilonijOverKb,
    const List<List<scalar>>& deltaij,
    const List<List<scalar>>& Mij,
    const List<List<scalar>>& sigmaij
)
{
    Ymd_             = Ymd;
    Xmd_             = Xmd;
    epsilonijOverKb_ = epsilonijOverKb;
    deltaij_         = deltaij;
    Mij_             = Mij;
    sigmaij_         = sigmaij;
}


// functions to access coefficients
template<class Thermo>
inline Foam::scalar Foam::sutherlandKineticTransport<Thermo>::epsilonOverKb() const
{
    return epsilonOverKb_;
}

template<class Thermo>
inline Foam::scalar Foam::sutherlandKineticTransport<Thermo>::sigma() const
{
    return sigma_;
}

template<class Thermo>
inline Foam::scalar Foam::sutherlandKineticTransport<Thermo>::miui() const
{
    return miui_;
}

template<class Thermo>
inline Foam::scalar Foam::sutherlandKineticTransport<Thermo>::polar() const
{
    return polar_;
}

template<class Thermo>
inline Foam::scalar Foam::sutherlandKineticTransport<Thermo>::Kb() const
{
    return Kb_;
}


template<class Thermo>
inline Foam::List<Foam::scalar>
Foam::sutherlandKineticTransport<Thermo>::Ymd() const
{
    return Ymd_;
}

template<class Thermo>
inline Foam::List<Foam::scalar>
Foam::sutherlandKineticTransport<Thermo>::Xmd() const
{
    return Xmd_;
}

template<class Thermo>
inline Foam::List<Foam::List<Foam::scalar>>
Foam::sutherlandKineticTransport<Thermo>::epsilonijOverKb() const
{
    return epsilonijOverKb_;
}

template<class Thermo>
inline Foam::List<Foam::List<Foam::scalar>>
Foam::sutherlandKineticTransport<Thermo>::deltaij() const
{
    return deltaij_;
}

template<class Thermo>
inline Foam::List<Foam::List<Foam::scalar>>
Foam::sutherlandKineticTransport<Thermo>::Mij() const
{
    return Mij_;
}

template<class Thermo>
inline Foam::List<Foam::List<Foam::scalar>>
Foam::sutherlandKineticTransport<Thermo>::sigmaij() const
{
    return sigmaij_;
}

// ************************************************************************* //
