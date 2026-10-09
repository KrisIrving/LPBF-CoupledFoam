// Optional spatial observations for a uniform Cartesian synthetic test mesh.
if (runTime.controlDict().lookupOrDefault<Switch>("spatialDiagnostics", false))
{
    const scalar condensedVolume = gSum(condensate.primitiveField()*mesh.V().field());
    const scalar condensedResistivity = gSum
    (
        electrical_resistivity.primitiveField()*condensate.primitiveField()*mesh.V().field()
    )/max(condensedVolume, VSMALL);
    if (Pstream::master())
    {
        const int oldPrecision = Info().precision();
        Info().precision(15);
        Info<< "M1_OPTICS time=" << runTime.value()
            << " condensedWeightedResistivityOhmM=" << condensedResistivity << nl;
        Info().precision(oldPrecision);
    }
    scalarField alphaSum(mesh.nCells(), Zero);
    for (const phaseModel& phase : mixture.phases())
    {
        alphaSum += phase.primitiveField();
        const scalar volume = gSum(phase.primitiveField()*mesh.V().field());
        scalar momentX = 0;
        forAll(phase, celli)
        {
            momentX += phase[celli]*mesh.C()[celli].x()*mesh.V()[celli];
        }
        reduce(momentX, sumOp<scalar>());
        const scalar minimum = gMin(phase.primitiveField());
        const scalar maximum = gMax(phase.primitiveField());
        if (Pstream::master())
        {
            const int oldPrecision = Info().precision();
            Info().precision(15);
            Info<< "M1_SPATIAL time=" << runTime.value()
                << " phase=" << phase.name()
                << " alphaMin=" << minimum << " alphaMax=" << maximum
                << " alphaVolumeM3=" << volume
                << " alphaMomentXM4=" << momentX << nl;
            Info().precision(oldPrecision);
        }
    }
    const scalar sumMin = gMin(alphaSum);
    const scalar sumMax = gMax(alphaSum);
    scalar analyticError = -1;
    if (runTime.controlDict().found("spatialReference"))
    {
        const dictionary& reference = runTime.controlDict().subDict("spatialReference");
        const scalar dx = reference.get<scalar>("cellWidthX");
        const scalar shift = reference.get<scalar>("velocityX")*runTime.value();
        const scalar lo = reference.get<scalar>("slabMinX")+shift;
        const scalar hi = reference.get<scalar>("slabMaxX")+shift;
        if (dx <= 0) FatalErrorInFunction << "cellWidthX must be positive" << exit(FatalError);
        scalar errorVolume = 0;
        for (const phaseModel& phase : mixture.phases())
        {
            if (phase.name() != "metal1") continue;
            forAll(phase, celli)
            {
                const scalar x = mesh.C()[celli].x();
                const scalar exact = max(scalar(0), min(scalar(1),
                    (min(x+0.5*dx, hi)-max(x-0.5*dx, lo))/dx));
                errorVolume += mag(phase[celli]-exact)*mesh.V()[celli];
            }
        }
        reduce(errorVolume, sumOp<scalar>());
        analyticError = errorVolume/gSum(mesh.V().field());
    }
    if (Pstream::master())
    {
        const int oldPrecision = Info().precision();
        Info().precision(15);
        Info<< "M1_INTERFACE time=" << runTime.value()
            << " alphaSumMin=" << sumMin << " alphaSumMax=" << sumMax
            << " analyticMeanAbsError=" << analyticError << nl;
        Info().precision(oldPrecision);
    }
}
