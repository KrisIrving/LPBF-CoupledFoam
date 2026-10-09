tmp<volScalarField> rAU;
if (restoringContinuumState)
{
    Info<< "M1_RESTART_STATE time=" << runTime.timeName()
        << " preservedFields=phi,rho,epsilon1 skippedInitialZeroDivergenceProjection=true"
        << " epsilonMin=" << gMin(epsilon1.primitiveField())
        << " epsilonMax=" << gMax(epsilon1.primitiveField())
        << endl;
}

if (correctPhi)
{
    rAU = new volScalarField
    (
        IOobject
        (
            "rAU",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rAU", dimTime/dimDensity, 1)
    );

    if (!restoringContinuumState)
    {
        #include "correctPhi.esi.H"
    }
}
else if (!restoringContinuumState)
{
    CorrectPhi
    (
        U,
        phi,
        p_rgh,
        dimensionedScalar("rAUf", dimTime/rho.dimensions(), 1),
        geometricZeroField(),//deallocation error if mDot here?!?!
        pimple
    );

// Info<<"! Corrected PHI for MDot !"<<endl;

    #include "continuityErrs.H"
}
