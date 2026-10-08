// Optional M1 observations. No equations or source terms are modified.
// Include with diagnosticSolver and diagnosticMassFlux in scope.
if (runTime.controlDict().lookupOrDefault<Switch>("continuumDiagnostics", false))
{
    const scalarField& diagnosticVolumes = mesh.V().field();
    const scalar diagnosticMass =
        gSum(rho.primitiveField()*diagnosticVolumes);
    const scalar diagnosticKineticEnergy =
        0.5*gSum(rho.primitiveField()*magSqr(U.primitiveField())*diagnosticVolumes);
    const scalar diagnosticLaserPower =
        gSum(laser.deposition().primitiveField()*diagnosticVolumes);
    const scalar diagnosticTmin = gMin(T.primitiveField());
    const scalar diagnosticTmax = gMax(T.primitiveField());
    const scalar diagnosticMaxSpeed = gMax(mag(U.primitiveField()));

    scalar diagnosticOutwardMassFlux = 0.0;
    forAll(diagnosticMassFlux.boundaryField(), patchi)
    {
        // Processor and periodic interfaces are internal to the global domain.
        if (!diagnosticMassFlux.boundaryField()[patchi].coupled())
        {
            diagnosticOutwardMassFlux += sum(diagnosticMassFlux.boundaryField()[patchi]);
        }
    }
    reduce(diagnosticOutwardMassFlux, sumOp<scalar>());

    if (Pstream::master())
    {
        const int diagnosticOldPrecision = Info().precision();
        Info().precision(12);
        Info<< "M1_CONTINUUM solver=" << diagnosticSolver
            << " time=" << runTime.value()
            << " dt=" << runTime.deltaTValue()
            << " massKg=" << diagnosticMass
            << " outwardMassFluxKgPerS=" << diagnosticOutwardMassFlux
            << " kineticEnergyJ=" << diagnosticKineticEnergy
            << " absorbedLaserPowerW=" << diagnosticLaserPower
            << " TminK=" << diagnosticTmin
            << " TmaxK=" << diagnosticTmax
            << " maxSpeedMS=" << diagnosticMaxSpeed
            << nl;
        Info().precision(diagnosticOldPrecision);
    }
}
