// Phase inventories include EOS density changes; they are not raw source rates.
if (runTime.controlDict().lookupOrDefault<Switch>("continuumDiagnostics", false))
{
    for (const phaseModel& diagnosticPhase : mixture.phases())
    {
        const tmp<volScalarField> diagnosticPhaseRho = diagnosticPhase.thermo().rho();
        const scalar diagnosticPhaseMass = gSum
        (
            diagnosticPhase.primitiveField()
           *diagnosticPhaseRho().primitiveField()*mesh.V().field()
        );
        const scalar diagnosticAlphaVolume =
            gSum(diagnosticPhase.primitiveField()*mesh.V().field());
        if (Pstream::master())
        {
            const int diagnosticPrecision = Info().precision();
            Info().precision(12);
            Info<< "M1_PHASE time=" << runTime.value()
                << " phase=" << diagnosticPhase.name()
                << " massKg=" << diagnosticPhaseMass
                << " alphaVolumeM3=" << diagnosticAlphaVolume << nl;
            Info().precision(diagnosticPrecision);
        }
    }
}
