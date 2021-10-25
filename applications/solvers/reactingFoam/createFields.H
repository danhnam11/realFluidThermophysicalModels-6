#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;
autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
     Info << "Nam - Check create pThermo" << endl; //Nam 
psiReactionThermo& thermo = pThermo();
     Info << "Nam - Check creat thermo = pThermo" << endl; //Nam 
thermo.validate(args.executable(), "h", "e");

     Info << "Nam - Check mixture is called before making composition" << endl; //Nam 

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

     Info << "Nam - Check mixture is called after making composition" << endl; //Nam 

const word inertSpecie(thermo.lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species() << exit(FatalIOError);
}

     Info << "Nam - Create rho in solver" << endl; //Nam 

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ, // Nam April 6th 2020
        IOobject::AUTO_WRITE // Nam April 6th 2020
    ),
    thermo.rho()
);

// Generate new field to validate Transport model - Nam April 6th 2020 
// generate dynamic viscosity field "mu" [kg/ms] 

     Info << "Nam - Create mu in solver" << endl;  //Nam

volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.mu()
);
// generate thermal conductivity field"kappa" [W/mK] 
volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.kappa()
);

// generate idea gas heat capacity at constant volume "Cv" [J/kgK] 
volScalarField Cv
(
    IOobject
    (
        "Cv",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.Cv()
);

// ------------------------------------------------------------- 

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField& p = thermo.p();

#include "compressibleCreatePhi.H"

pressureControl pressureControl(p, rho, pimple.dict(), false);

mesh.setFluxRequired(p.name());

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating reaction model\n" << endl;
autoPtr<CombustionModel<psiReactionThermo>> reaction
(
    CombustionModel<psiReactionThermo>::New(thermo, turbulence())
);

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

volScalarField Qdot
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Qdot", dimEnergy/dimVolume/dimTime, 0.0)
);

#include "createMRF.H"
#include "createFvOptions.H"
