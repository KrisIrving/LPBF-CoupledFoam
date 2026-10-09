{
    const Switch condensedPhaseOptics
    (
        runTime.controlDict().lookupOrDefault<Switch>("condensedPhaseOptics", false)
    );
    alphasum *= 0.0;
    condensate = 1.0;

    // refine_test*= 0.0;

    TSolidus *= 0.0;
    TLiquidus *= 0.0;
    LatentHeat *= 0.0;

    // electrical resistivity of condensed mixture/
    electrical_resistivity *= 0.0;
    eleccond *= 0.0;
    magperm *= 0.0;

    forAllConstIter(PtrDictionary<phaseModel>, mixture.phases(), iter)
    {
        const phaseModel& phase1 = iter();
        const volScalarField& alpha1 = phase1;

        IOdictionary phasedictionary
        (
            IOobject
            (
                "thermophysicalProperties." + phase1.name(),
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED
            )
        );

        dimensionedScalar phaseTSolidus
        (
            "phaseTSolidus",
            dimensionSet(0, 0, 0, 1, 0),
            readScalar(phasedictionary.lookup("TSolidus"))
        );

        dimensionedScalar phaseTLiquidus
        (
            "phaseTLiquidus",
            dimensionSet(0, 0, 0, 1, 0),
            readScalar(phasedictionary.lookup("TLiquidus"))
        );

        dimensionedScalar phaseLH
        (
            "phaseLH",
            dimensionSet(0, 2, -2, 0, 0),
            readScalar(phasedictionary.lookup("LatentHeat"))
        );

        dimensionedScalar phaseconductivity
        (
            "phaseconductivity",
            dimensionSet(-1, -3, 3, 0, 0),
            readScalar(phasedictionary.lookup("elec_conductivity"))
        );

        dimensionedScalar phasepermeability
        (
            "phasepermeability",
            dimensionSet(1, 1, -2, 0, 0),
            readScalar(phasedictionary.lookup("magnetic_permeability"))
        );


        TSolidus += alpha1*phaseTSolidus;
        TLiquidus += alpha1*phaseTLiquidus;
        LatentHeat += alpha1*phaseLH;
        eleccond += alpha1*phaseconductivity;
        magperm += alpha1*phasepermeability;


        std::string str2("vapour");

        std::string phasename(phase1.name());

        const unsigned long res = phasename.find(str2);

        if (res != string::npos || phase1.name() == "air")
        {
            condensate -= alpha1;
            if (!condensedPhaseOptics)
            {
                electrical_resistivity += alpha1/phaseconductivity;
            }
        }
        else if (condensedPhaseOptics)
        {
            if (phaseconductivity.value() <= 0)
            {
                FatalErrorInFunction << "Condensed-phase conductivity must be positive"
                    << exit(FatalError);
            }
            electrical_resistivity += max(alpha1, scalar(0))/phaseconductivity;
        }

        alphasum += alpha1;
    }

    condensate = max(min(condensate, 1.0), 0.0);
    if (condensedPhaseOptics)
    {
        // Resistivity of the condensed mixture, not diluted by gas fraction.
        // Single-metal limit: eta=1/sigma. Default preserves upstream behavior.
        electrical_resistivity /= max(condensate, scalar(1e-12));
    }


    electrical_resistivity.correctBoundaryConditions();
    magperm.correctBoundaryConditions();
    eleccond.correctBoundaryConditions();
    LatentHeat.correctBoundaryConditions();
    TLiquidus.correctBoundaryConditions();
    TSolidus.correctBoundaryConditions();


    // Update the cell size fields
    const faceList& ff = mesh.faces();
    const pointField& pp = mesh.points();

    const vectorField& CI = mesh.C();
    scalarField& yDimI = laser.yDim();

    forAll(CI, celli)
    {
        // vector XYZ = mesh.C()[celli];
        // xcoord[celli]=XYZ.x();
        // zcoord[celli]=XYZ.z();

        const cell& cc = mesh.cells()[celli];
        labelList pLabels(cc.labels(ff));
        pointField pLocal(pLabels.size(), vector::zero);

        forAll(pLabels, pointi)
        {
            pLocal[pointi] = pp[pLabels[pointi]];
        }

        // xDim[celli] =
        //     Foam::max(pLocal & vector(1,0,0))
        //   - Foam::min(pLocal & vector(1,0,0));
        yDimI[celli] =
            Foam::max(pLocal & vector(0,1,0))
          - Foam::min(pLocal & vector(0,1,0));
        // zDim[celli] =
        //     Foam::max(pLocal & vector(0,0,1))
        //   - Foam::min(pLocal & vector(0,0,1));
    }

    laser.yDim().correctBoundaryConditions();

    metalaverage = fvc::average(condensate);

    // laser.refineFlag() *= 0.0;
    // laser.refineFlag() += metalaverage;

    const volVectorField gradcondensate(fvc::grad(condensate));

    n_filtered = fvc::grad(condensate)/((mag(gradcondensate) + deltaN));
}
