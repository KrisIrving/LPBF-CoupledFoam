// End-of-step pressure response. Sources are logged separately at alpha solves.
if (runTime.controlDict().lookupOrDefault<Switch>("interfaceDiagnostics", false))
{
    const scalar pressureMean = gSum(p.primitiveField()*mesh.V().field())/gSum(mesh.V().field());
    const scalar pressureMin = gMin(p.primitiveField());
    const scalar pressureMax = gMax(p.primitiveField());
    if (Pstream::master())
    {
        const int savedPrecision = Info().precision();
        Info().precision(15);
        Info<< "M1C_FLOW time=" << runTime.value() << " pressureMeanPa=" << pressureMean
            << " pressureMinPa=" << pressureMin << " pressureMaxPa=" << pressureMax << nl;
        Info().precision(savedPrecision);
    }
}
