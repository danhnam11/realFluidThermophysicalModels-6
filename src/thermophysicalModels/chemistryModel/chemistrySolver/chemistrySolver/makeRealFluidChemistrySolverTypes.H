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

    Macro to make chemistry solvers for real-fluid thermophysicalModels library

    by:      
    Nam Danh Nguyen - PhD.st
    Advisor: Prof. Chun Sang Yoo 
    Combustion & Propulsion Lab - Dept. of Mech. Engineering
    Ulsan Institute of Science and Technology (UNIST) - Korea	

\*---------------------------------------------------------------------------*/

#ifndef makeRealFluidChemistrySolverTypes_H
#define makeRealFluidChemistrySolverTypes_H

#include "chemistrySolver.H"
#include "SRKchungTakaStandardChemistryModel.H"
#include "SRKchungKineticStandardChemistryModel.H"
#include "idKineticStandardChemistryModel.H"
#include "PRchungTakaStandardChemistryModel.H"
#include "PRchungKineticStandardChemistryModel.H"
//#include "TDACChemistryModel.H"

#include "noChemistrySolver.H"
#include "EulerImplicit.H"
#include "ode.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSRKchungTakaChemistrySolverType(SS, Comp, Thermo)                  \
                                                                               \
    typedef SS<SRKchungTakaStandardChemistryModel<Comp, Thermo>> SS##Comp##Thermo; \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SS##Comp##Thermo,                                                      \
        (#SS"<" + word(SRKchungTakaStandardChemistryModel                      \
                <Comp, Thermo>::typeName_()) + "<" + word(Comp::typeName_())   \
                + "," + Thermo::typeName() + ">>").c_str(),                    \
        0                                                                      \
    );                                                                         \
                                                                               \
    BasicChemistryModel<Comp>::                                                \
        add##thermo##ConstructorToTable<SS##Comp##Thermo>                      \
        add##SS##Comp##Thermo##thermo##ConstructorTo##BasicChemistryModel##Comp\
##Table_; 
                                                                               
#define makeSRKchungTakaChemistrySolverTypes(Comp, Thermo)                     \
                                                                               \
    makeSRKchungTakaChemistrySolverType                                        \
    (                                                                          \
        noChemistrySolver,                                                     \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeSRKchungTakaChemistrySolverType                                        \
    (                                                                          \
        EulerImplicit,                                                         \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeSRKchungTakaChemistrySolverType                                        \
    (                                                                          \
        ode,                                                                   \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \



#define makeSRKchungKineticChemistrySolverType(SS, Comp, Thermo)               \
                                                                               \
    typedef SS<SRKchungKineticStandardChemistryModel<Comp, Thermo>> SS##Comp##Thermo; \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SS##Comp##Thermo,                                                      \
        (#SS"<" + word(SRKchungKineticStandardChemistryModel                   \
                <Comp, Thermo>::typeName_()) + "<" + word(Comp::typeName_())   \
                + "," + Thermo::typeName() + ">>").c_str(),                    \
        0                                                                      \
    );                                                                         \
                                                                               \
    BasicChemistryModel<Comp>::                                                \
        add##thermo##ConstructorToTable<SS##Comp##Thermo>                      \
        add##SS##Comp##Thermo##thermo##ConstructorTo##BasicChemistryModel##Comp\
##Table_; 
                                                                               
#define makeSRKchungKineticChemistrySolverTypes(Comp, Thermo)                  \
                                                                               \
    makeSRKchungKineticChemistrySolverType                                     \
    (                                                                          \
        noChemistrySolver,                                                     \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeSRKchungKineticChemistrySolverType                                     \
    (                                                                          \
        EulerImplicit,                                                         \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeSRKchungKineticChemistrySolverType                                     \
    (                                                                          \
        ode,                                                                   \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \



#define makePRchungTakaChemistrySolverType(SS, Comp, Thermo)                   \
                                                                               \
    typedef SS<PRchungTakaStandardChemistryModel<Comp, Thermo>> SS##Comp##Thermo; \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SS##Comp##Thermo,                                                      \
        (#SS"<" + word(PRchungTakaStandardChemistryModel                       \
                <Comp, Thermo>::typeName_()) + "<" + word(Comp::typeName_())   \
                + "," + Thermo::typeName() + ">>").c_str(),                    \
        0                                                                      \
    );                                                                         \
                                                                               \
    BasicChemistryModel<Comp>::                                                \
        add##thermo##ConstructorToTable<SS##Comp##Thermo>                      \
        add##SS##Comp##Thermo##thermo##ConstructorTo##BasicChemistryModel##Comp\
##Table_; 
                                                                               
#define makePRchungTakaChemistrySolverTypes(Comp, Thermo)                      \
                                                                               \
    makePRchungTakaChemistrySolverType                                         \
    (                                                                          \
        noChemistrySolver,                                                     \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makePRchungTakaChemistrySolverType                                         \
    (                                                                          \
        EulerImplicit,                                                         \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makePRchungTakaChemistrySolverType                                         \
    (                                                                          \
        ode,                                                                   \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \


#define makePRchungKineticChemistrySolverType(SS, Comp, Thermo)                \
                                                                               \
    typedef SS<PRchungKineticStandardChemistryModel<Comp, Thermo>> SS##Comp##Thermo; \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SS##Comp##Thermo,                                                      \
        (#SS"<" + word(PRchungKineticStandardChemistryModel                    \
                <Comp, Thermo>::typeName_()) + "<" + word(Comp::typeName_())   \
                + "," + Thermo::typeName() + ">>").c_str(),                    \
        0                                                                      \
    );                                                                         \
                                                                               \
    BasicChemistryModel<Comp>::                                                \
        add##thermo##ConstructorToTable<SS##Comp##Thermo>                      \
        add##SS##Comp##Thermo##thermo##ConstructorTo##BasicChemistryModel##Comp\
##Table_; 
                                                                               
#define makePRchungKineticChemistrySolverTypes(Comp, Thermo)                   \
                                                                               \
    makePRchungKineticChemistrySolverType                                      \
    (                                                                          \
        noChemistrySolver,                                                     \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makePRchungKineticChemistrySolverType                                      \
    (                                                                          \
        EulerImplicit,                                                         \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makePRchungKineticChemistrySolverType                                      \
    (                                                                          \
        ode,                                                                   \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \



#define makeidKineticChemistrySolverType(SS, Comp, Thermo)                     \
                                                                               \
    typedef SS<idKineticStandardChemistryModel<Comp, Thermo>> SS##Comp##Thermo;\
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SS##Comp##Thermo,                                                      \
        (#SS"<" + word(idKineticStandardChemistryModel                         \
                <Comp, Thermo>::typeName_()) + "<" + word(Comp::typeName_())   \
                + "," + Thermo::typeName() + ">>").c_str(),                    \
        0                                                                      \
    );                                                                         \
                                                                               \
    BasicChemistryModel<Comp>::                                                \
        add##thermo##ConstructorToTable<SS##Comp##Thermo>                      \
        add##SS##Comp##Thermo##thermo##ConstructorTo##BasicChemistryModel##Comp\
##Table_; 
                                                                               
#define makeidKineticChemistrySolverTypes(Comp, Thermo)                        \
                                                                               \
    makeidKineticChemistrySolverType                                           \
    (                                                                          \
        noChemistrySolver,                                                     \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeidKineticChemistrySolverType                                           \
    (                                                                          \
        EulerImplicit,                                                         \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeidKineticChemistrySolverType                                           \
    (                                                                          \
        ode,                                                                   \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
