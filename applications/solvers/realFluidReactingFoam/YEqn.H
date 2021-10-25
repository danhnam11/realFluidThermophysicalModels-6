
tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

//Info << "\ninert index = " << inertIndex << " and Name = " << speciesName[inertIndex] << endl;
{
    //- Reset diffusion velocity to zero for each time step
    forAll(YVi, i)
    {
        YVi[i] *= 0.;
    }
    //---
    reaction->correct();
    Qdot = reaction->Qdot();
    volScalarField Yt(0.0*Y[0]);
    forAll(Y, i)
    {
        // calculate correction term
        sumDiffusionCorrections1 *= 0.;
        sumDiffusionCorrections2 *= 0.; 
        forAll(Y, k)
        {
            sumDiffusionCorrections1 += fvc::laplacian(Y[i]*rho*Dimix[k], Y[k]);
            sumDiffusionCorrections2 += fvc::laplacian(Y[i]*rho*Dimix[k]*Y[k]/Wmix, Wmix);
        }

        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(rho*Dimix[i], Yi)        
             == 
                reaction->R(Yi)
              + fvc::laplacian(rho*Dimix[i]*Yi/Wmix, Wmix)   
              - sumDiffusionCorrections1                     
              - sumDiffusionCorrections2                     
              + fvOptions(rho, Yi)
            );
 
            YiEqn.relax();
            fvOptions.constrain(YiEqn); 
            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);
            Yi.max(0.0);
            Yt += Yi;
        }
    }
    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);

    //-Calculate diffusion velocity of inert species
    volVectorField YVt(0.0*YVi[0]); 
    forAll(YVi, i)
    {
        //- Update diffusion velocity for all species
        YVi[i] = (-(Dimix[i]/Wmix)*Wmix*fvc::grad(Y[i])-(Dimix[i]/Wmix)*Y[i]*fvc::grad(Wmix));
        if (i != inertIndex && composition.active(i))
        {
            YVt += YVi[i];
        }
    }
    YVi[inertIndex] = -YVt; 
}
