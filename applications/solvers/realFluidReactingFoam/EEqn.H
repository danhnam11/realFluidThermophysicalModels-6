{
    volScalarField& he = thermo.he();

    //-Reset and then calculate Heat by diffusion terms-for each time step
    sumHeatDiffusion *= 0.;
    sumHeatDiffusion2 *= 0.;

    forAll(Y, k)
    {
        sumHeatDiffusion += fvc::laplacian((kappa/Cp)*hei[k], Y[k]); 
        sumHeatDiffusion2 += fvc::div(hei[k]*rho*YVi[k]); 
    }
    // -------

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)  
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(kappa/Cp, he) 
     ==
      Qdot
      - sumHeatDiffusion     
      - sumHeatDiffusion2    
      + fvOptions(rho, he)
    );
    
    EEqn.relax();
    fvOptions.constrain(EEqn);
    EEqn.solve();
    fvOptions.correct(he);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}


