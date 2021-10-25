# A real-fluid based thermophysicalModels library OF-6
## General Information
An updated _thermophysicalModels_ libary of OpenFOAM 6.0. with real-fluid models for reacting flow simulations at high pressure. Readers are referred to documentations provided in _documentations_ directory for the detail implementation and extension guide. They are written for the development in OpenFOAM-6 but they can be referred for the development in other versions.
## List of implemented real-fluid models in the new libary
- Modified Soave-Redlich-Kwong (SRK) model for equation of state [1, 2].
- Peng-Robinson (PR) model for equation of state [3]
- JANAF-based model for real-fluid thermodynamic properties.
- Chung's model (1988) for dynamic viscosity and thermal conductivity [4].
- Mixture averaged model for mass diffusivity of individual species in a mixture in which the binary diffusion coefficients are based on Fuller's model and Takahashi correction at high pressure [5].
- Mixture averaged model for mass diffusivity of individual species in a mixture in which the binary diffusion coefficients are based on Standard Kinetic Theory [6].
## Combinations of implemented models available in the library
The runtime names of thermophysical models need to be specified correctly in the _thermoType_ dictionary and _chemistryType_ dictionary based on the combinations of implemented models available in the libary. For example:

In the _constant/thermophysicalProperties_ file: 

	thermoType
	{
		type            hePsiThermo;                   //(1)
		mixture         SRKchungTakaReactingMixture;   //(2)
		transport       chungTaka;                     //(3)
		thermo          rfJanaf;                       //(4)
		energy          sensibleEnthalpy;              //(7)
		equationOfState soaveRedlichKwong;             //(5)
		specie          rfSpecie;                      //(6)
	}

In the _constant/chemistryProperties_ file: 

	chemistryType
	{
		solver   ode;                                  // either ode/EulerImplicit/none;
		method   SRKchungTakaStandard;                 //(8)
	}

There are 20 options (combinations) are available for reacting flow simulations in the library as follows:

| No | type(1)| mixture(2) | transport(3) | thermo(4) | EoS(5) | specie(6) | Energy(7) | method(8)      |
| :- | :----- |:-----------| :----------- | :-------- | :----- | :-------- | :-------- | :------------- |
|1   | Psi/Rho| SRK-C-Ta   | C-Ta         | rfJ       | SRK    | rfSp      | sens/Int  | SRK-C-Ta-Stand |
|2   | Psi/Rho| SRK-C-Ki   | C-Ki         | rfJ       | SRK    | rfSp      | sens/Int  | SRK-C-Ki-Stand |
|3   | Psi/Rho| PR-C-Ta    | C-Ta         | rfJ       | PR     | rfSp      | sens/Int  | PR-C-Ta-Stand  |
|4   | Psi/Rho| PR-C-Ki    | C-Ki         | rfJ       | PR     | rfSp      | sens/Int  | PR-C-Ki-Stand  |
|5   | Psi/Rho| id-Ki      | su-Ki        | J         | per    | sp        | sens/Int  | id-Ki-Stand    |

where abbreviations of model names are:

(1) Psi: _hePsiThermo_; Rho: _heRhoThermo_.
(2) SRK-C-Ta: _SRKchungTakaReactingMixture_; SRK-C-Ki: _SRKchungKineticReactingMixture_; PR-C-Ta: _PRchungTakaReactingMixture_; PR-C-Ki: _PRchungKineticReactingMixture_; id-Ki: _idKineticReactingMixture_.
(3) C-Ta: _chungTaka_; C-Ki: _chungKinetic_; su-Ki: _sutherlandKinetic_.
(4) rfJ: _rfJanaf_; J: _janaf_.
(5) SRK: _soaveRedlichKwong_; PR: _Peng_; per: _perfectGas_.
(6) rfSp: _rfSpecie_; sp: _specie_.
(7) sens: _sensibleEnthalpy_; Int: _sensibleInternalEnergy_.
(8) SRK-C-Ta-Stand: _SRKchungTakaStandard_; SRK-C-Ki-Stand: _SRKchungKineticStandard_; PR-C-Ta-Stand: _PRchungTakaStandard_; PR-C-Ki-Stand: _PRchungKineticStandard_; id-Ki-Stand: _idKineticStandard_.


## Installation
- Since the new library is developed based on OpenFOAM 6.0 in Linux operating systems, the complete installation of OpenFOAM 6.0 framework is required. 
- Prepare a directory on your system, e.g., _yourDirectory_:

		mkdir ~/OpenFOAM/yourDirectory/
		cd ~/OpenFOAM/yourDirectory/	
- Download source files using git: 

		git clone https://github.com/danhnam11/realFluidThermophysicalModels-6.git

- Specify the path of your _src_ directory to an environment variable, _LIB_REALFLUID_SRC_. For example:

		echo "export LIB_REALFLUID_SRC=~/OpenFOAM/yourDirectory/realFluidThermophysicalModels-6/src/" >> ~/.bashrc
		source ~/.bashrc
- To compile the library and applications, go to _realFluidThermophysicalModels_ directory and run the _Allwmake_ script:

		cd ~/OpenFOAM/yourDirectory/realFluidThermophysicalModels-6/
		./Allwmake

- Now the real-fluid based _thermophysicalModels_ library and _realFluidReactingFoam_ solver are ready to use. The library and solver are stored at _$FOAM_USER_LIBBIN_ and _$FOAM_USER_APPBIN_ directory.

## Using real-fluid models in a solver 
- See _implementationGuide.pdf_ file in documentations.
- Applications of interest that adopt _thermophysicalModels_, _combustionModels_, _compressibleTurbulenceModels_ libraries should be recompiled with a proper path linking to the new ones to avoid a _segmentation fault_ error since new created object files of those libraries have been replaced to the original ones in OpenFOAM.
- A function returing mass diffusivity has not been implemented yet for rho-based system.
- Readers are referred to our paper for the validation of the new library.

## Tutorials
Tutorials for generating thermodynamic and transport properties of a pure species, mixture and for laminar non-premixed counterflow flames of CH4 versus O2/CO2 are available in _tutorial_ directory.

	cd ~/OpenFOAM/yourDirectory/realFluidThermophysicalModels-6/tutorials/

## Notes:
Should you find bugs or have suggestions on how to make the code better, please post on cfd-online using following thread: http://www.cfd-online.com/Forums/ (This will be update later after the paper is accepted). 


## Authors 
This package was developed at Combustion & Propulsion Lab., Dept. of Mech. Engineering, Ulsan National Institute of Science and Technology (UNIST), Korea (Prof. C.S. Yoo: https://csyoo.unist.ac.kr/). If you publish results that are obtained using this package, please cite our paper as follows:
- D. N. Nguyen, K. S. Jung, J. W. Shim, C. S. Yoo, Real-fluid thermophysicalModels: An OpenFOAM-based library for reacting flow simulations at high pressure, Computer Physics Communications (2021)(submitted).
## Reference
- [1] G. Soave, Equilibrium constants from a modified Redlich-Kwong equation of state, Chem. Eng. Sci. 27 (1972) 1197-1203.
- [2] D. Peng, D. Robinson, New two-equation of state, Ind. Eng. Chem. Fundam. 15(1976) 59-64. 
- [3] M. S. Graboski, T. E. Daubert, A modified Soave equation of state for phase equilibrium calculations. 1. Hydrocarbon systems, Ind. Eng. Chem. Process. Des. Dev. 17 (1978) 443-448.
- [4] T. C. Horng, M. Ajlan, L. L. Lee, K. E. Starling, M. Ajlan, Generalized multiparameter correlation for nonpolar and polar fluid transport properties, Ind. Eng. Chem. Res. 27 (1988) 671-679.
- [5] S. Takahashi, S. Takahashi, Preparation of a generalized chart for the diffusion coefficients of gases at high pressures, J. Chem. Eng. Japan 7 (1975) 417-420. 
- [6] R. J. Kee, F. M. Rupley, E. Meeks, J. A. Miller, CHEMKIN-III: a fortran chemical kinetics package for the analysis of gas-phase chemical and plasma kinetics, SAND96-8216 (1996). 
