// Observe continuity density versus EOS density without modifying either.
if (runTime.controlDict().lookupOrDefault<Switch>("couplingDiagnostics", false))
{
    const scalar continuityMass = gSum(rho.primitiveField()*mesh.V().field());
    const tmp<volScalarField> eosRho = mixture.rho();
    const scalar eosMass = gSum(eosRho().primitiveField()*mesh.V().field());
    const scalar domainVolume = gSum(mesh.V().field());
    const scalar pressureMean = gSum(p.primitiveField()*mesh.V().field())/domainVolume;
    const scalar volumeSource = gSum(vDot.primitiveField()*mesh.V().field());
    const scalar temperatureSource = gSum(mass_dot.primitiveField()*mesh.V().field());
    if (Pstream::master())
    {
        const int oldPrecision = Info().precision();
        Info().precision(15);
        Info<< "M1_COUPLING time=" << runTime.value()
            << " stage=" << couplingStage
            << " continuityMassKg=" << continuityMass
            << " eosMassKg=" << eosMass
            << " meanPressurePa=" << pressureMean
            << " volumeSourceM3PerS=" << volumeSource
            << " temperatureSourceKgKPerS=" << temperatureSource << nl;
        Info().precision(oldPrecision);
    }
}
