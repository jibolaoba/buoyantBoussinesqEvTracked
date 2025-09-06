/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

Application
    buoyantBoussinesqEvTracked

Group
    grpHeatTransferSolvers

Description
    Transient solver modified by Jibola Owolabi for buoyant turbulent flow with simplified solvent-solute
    model for RH-driven evaporation of aerosol particles in incompressible
    framework, with optional mesh motion and topology changes. Tracks droplet
    size reduction from initial to equilibrium diameter.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "radiationModel.H"
#include "CorrectPhi.H"
#include "basicKinematicCloud.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "OFstream.H"

// **Simplified Evaporation Model Function with Dynamic Density Update**
void updateDropletDiameter(
    Foam::KinematicParcel<Foam::particle>& particle,
    Foam::scalar deltaTime,
    Foam::scalar C,
    Foam::scalar RH,
    Foam::scalar initialD,
    const Foam::IOdictionary& kinematicCloudProperties
)
{
    // Read composition and densities from kinematicCloudProperties
    const scalar phi_solute_initial = kinematicCloudProperties.subDict("evaporationProperties").get<scalar>("phiSoluteInitial");
    const scalar phi_salts = kinematicCloudProperties.subDict("evaporationProperties").get<scalar>("phiSalts");
    const scalar phi_proteins = kinematicCloudProperties.subDict("evaporationProperties").get<scalar>("phiProteins");
    const scalar densitySalts = kinematicCloudProperties.subDict("evaporationProperties").get<scalar>("densitySalts");
    const scalar densityProteins = kinematicCloudProperties.subDict("evaporationProperties").get<scalar>("densityProteins");
    const scalar densityWater = kinematicCloudProperties.subDict("evaporationProperties").get<scalar>("densityWater");
    const scalar kappa = kinematicCloudProperties.subDict("evaporationProperties").get<scalar>("kappa");
    const scalar phi_water = 1.0 - phi_solute_initial;

    // Validate composition
    if (mag(phi_solute_initial - (phi_salts + phi_proteins)) > 1e-6)
    {
        FatalErrorInFunction
            << "phi_solute_initial must equal phi_salts + phi_proteins"
            << abort(FatalError);
    }
    if (mag(phi_water + phi_solute_initial - 1.0) > 1e-6)
    {
        FatalErrorInFunction
            << "phi_water + phi_solute_initial must equal 1.0"
            << abort(FatalError);
    }

    // Calculate initial droplet density (mass-weighted)
    //scalar densityDroplet = phi_water * densityWater + phi_salts * densitySalts + phi_proteins * densityProteins;

    // Calculate initial mass and solute volumes
    //scalar initialMass = (M_PI / 6.0) * Foam::pow3(initialD) * densityDroplet;
   // scalar soluteVolume = (phi_salts * initialMass / densitySalts) + (phi_proteins * initialMass / densityProteins);
    
    //to simplify
    //scalar V_droplet = (M_PI / 6.0) * Foam::pow3(initialD);
    //scalar soluteVolume = phi_solute_initial * V_droplet;
    //scalar D_dry = Foam::cbrt(6.0 * soluteVolume / M_PI);

   //  Wet droplet density (unchanged)
      scalar densityDroplet = phi_water * densityWater + phi_salts * densitySalts + phi_proteins * densityProteins;

   //  Volume and mass
      scalar V_droplet = (M_PI/6.0) * pow3(initialD);
      scalar m_droplet = densityDroplet * V_droplet;

   //  Mass fractions
      scalar phi_salts_m = phi_salts * densitySalts / densityDroplet;
      scalar phi_prot_m = phi_proteins * densityProteins / densityDroplet;

   //  Solute mass & volume
      scalar m_salts  = phi_salts_m  * m_droplet;
      scalar m_prot   = phi_prot_m   * m_droplet;

      scalar soluteVolume = m_salts / densitySalts + m_prot / densityProteins;
 
   //  Dry diameter
      scalar D_dry = Foam::cbrt(6.0 * soluteVolume / M_PI);

    // Calculate equilibrium diameter
    scalar growthFactor = Foam::pow(1.0 + kappa * RH / (1.0 - RH), 1.0/3.0);
    scalar Deq = D_dry * growthFactor;

    // Current size and evaporation
    scalar currentD = particle.d();
    scalar V_current = (M_PI / 6.0) * Foam::pow3(currentD);
    scalar V_eq = (M_PI / 6.0) * Foam::pow3(Deq);
    scalar surfaceArea = M_PI * Foam::sqr(currentD);
    scalar mDot = C * (1.0 - RH) * surfaceArea;
    scalar deltaMass = mDot * deltaTime;
    scalar deltaV = deltaMass / densityWater;

    // Update volume and diameter
    scalar V_new = max(V_current - deltaV, V_eq);
    scalar newD = Foam::cbrt(6.0 * V_new / M_PI);
    particle.d() = newD;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for buoyant, turbulent flow with simplified solvent-solute "
        "model for RH-driven evaporation of aerosol particles."
    );

    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    #include "initContinuityErrs.H"
    #include "createClouds.H"

    turbulence->validate();

    // **Read Kinematic Cloud Properties**
    IOdictionary kinematicCloudProperties
    (
        IOobject
        (
            "kinematicCloudProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    // **Read Initial Diameters for Models**
    HashTable<dimensionedScalar> initialDiameters;
    const dictionary& injectionModels = kinematicCloudProperties.subDict("subModels").subDict("injectionModels");
    wordList modelNames = {"model03_1"};  // Default to single model
    forAll(modelNames, i)
    {
        const word& modelName = modelNames[i];
        dimensionedScalar initialD
        (
            "initialDiameter_" + modelName,
            dimLength,
            injectionModels.subDict(modelName)
                          .subDict("sizeDistribution")
                          .subDict("fixedValueDistribution")
                          .get<scalar>("value")
        );
        initialDiameters.insert(modelName, initialD);
    }
    
    dimensionedScalar evapConst
    (
        "evapConst",
        dimMass/(dimTime*dimArea),
        kinematicCloudProperties.subDict("evaporationProperties").get<scalar>("C")
    );

    scalar sampleInterval = kinematicCloudProperties.subDict("postProcessing").get<scalar>("sampleInterval");
    scalar nextWriteTime = 0.0;

    DynamicList<scalar> parcelInitialDiameters(0);

    autoPtr<OFstream> dropletSizeFile;
    if (Pstream::master())
    {
        mkDir(runTime.globalPath() / "postProcessing");
        dropletSizeFile.reset
        (
            new OFstream(runTime.globalPath() / "postProcessing" / "dropletSize.csv")
        );
        dropletSizeFile() << "Time(s),DInitialDiameter(um),Diameter(um)" << nl;
    }

    Info << "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info << "Time = " << runTime.timeName() << nl << endl;

        parcels.storeGlobalPositions();
        mesh.update();

        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.controlledUpdate();

                if (mesh.changing())
                {
                    MRF.update();
                    if (correctPhi)
                    {
                        phi = mesh.Sf() & Uf();
                        #include "correctPhi.H"
                        fvc::makeRelative(phi, U);
                    }
                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            parcels.evolve();

            if (parcels.nParcels() > parcelInitialDiameters.size())
            {
                parcelInitialDiameters.setSize(parcels.nParcels(), 0.0);
                label parcelIndex = 0;
                forAll(parcels, i)
                {
                    label modelIndex = (parcelIndex % modelNames.size());
                    word modelName = modelNames[modelIndex];
                    parcelInitialDiameters[i] = initialDiameters[modelName].value();
                    parcelIndex++;
                }
            }

            label parcelIdx = 0;
            for (auto& particle : parcels)
            {
                scalar initialD = parcelInitialDiameters[parcelIdx];
                updateDropletDiameter
                (
                    particle,
                    runTime.deltaT().value(),
                    evapConst.value(),
                    RH.value(),
                    initialD,
                    kinematicCloudProperties
                );
                parcelIdx++;
            }

            
         if (Pstream::master() && runTime.value() >= nextWriteTime)
            {
                // ✅ NEW: Added initial diameter logging per parcel
                label i = 0;
                for (auto& p : parcels)
                {
                    scalar initialD_um = parcelInitialDiameters[i] * 1e6;
                    scalar currentD_um = p.d() * 1e6;

                    dropletSizeFile() << runTime.value() << ","
                                      << initialD_um << ","
                                      << currentD_um << nl;

                    ++i;
                }
                dropletSizeFile().flush();
                nextWriteTime += sampleInterval;
            }



            #include "UEqn.H"
            #include "TEqn.H"

            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
