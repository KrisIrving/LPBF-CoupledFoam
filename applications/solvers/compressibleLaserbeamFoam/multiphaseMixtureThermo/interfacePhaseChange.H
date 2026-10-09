// Opt-in interface mass/energy contract. Legacy closure remains the default.
// Included before any phase is advanced, so both sources use identical states.
const dictionary& interfaceControls = mesh_.solverDict("alpha");
const word interfaceModel
(
    interfaceControls.lookupOrDefault<word>("interfacePhaseChangeModel", "legacy")
);
const bool interfacePhaseChange = interfaceModel != "legacy";
PtrList<volScalarField> interfaceAlphaSources
(
    interfacePhaseChange ? phases_.size() : 0
);
if (interfacePhaseChange)
{
    if
    (
        (interfaceModel != "prescribed" && interfaceModel != "kinetic")
     || !commonLatentHeatSource || !referenceInternalEnergySource
     || interfaceControls.get<label>("nAlphaSubCycles") != 1
    )
    {
        FatalErrorInFunction
            << "Interface candidate requires prescribed/kinetic, commonLatentHeatSource, "
            << "referenceInternalEnergySource and nAlphaSubCycles=1"
            << exit(FatalError);
    }
    const word liquidName(interfaceControls.get<word>("interfaceLiquid"));
    const word vapourName(interfaceControls.get<word>("interfaceVapour"));
    phaseModel* liquid = nullptr;
    phaseModel* vapour = nullptr;
    label liquidIndex = -1, vapourIndex = -1, index = 0;
    for (phaseModel& phase : phases_)
    {
        interfaceAlphaSources.set
        (
            index,
            new volScalarField
            (
                IOobject("interfaceSource_" + phase.name(), mesh_.time().timeName(), mesh_),
                mesh_, dimensionedScalar("zero", dimless/dimTime, 0)
            )
        );
        if (phase.name() == liquidName) { liquid = &phase; liquidIndex = index; }
        else if (phase.name() == vapourName) { vapour = &phase; vapourIndex = index; }
        else if (gMax(mag(phase.primitiveField())) > 1e-10)
        {
            FatalErrorInFunction << "Interface candidate currently supports a pure liquid/vapour pair; "
                << "active extra phase " << phase.name() << exit(FatalError);
        }
        ++index;
    }
    if (!liquid || !vapour || vapourName != liquidName + "vapour")
    {
        FatalErrorInFunction << "Missing/distinct liquid-vapour pair or unsupported naming"
            << exit(FatalError);
    }
    const scalar prescribedJ = interfaceControls.lookupOrDefault<scalar>("interfaceMassFlux", 0);
    const scalar accommodation = interfaceControls.lookupOrDefault<scalar>("interfaceAccommodation", 1);
    const scalar inventoryFraction = interfaceControls.lookupOrDefault<scalar>("interfaceInventoryFraction", 0.2);
    if (!std::isfinite(prescribedJ) || !std::isfinite(accommodation)
     || !std::isfinite(inventoryFraction) || accommodation < 0 || accommodation > 1
     || inventoryFraction <= 0 || inventoryFraction > 0.5)
    {
        FatalErrorInFunction << "Invalid accommodation or inventory fraction" << exit(FatalError);
    }
    const dictionary& offsets = phasedictionary.subDict("phaseEnergyOffsets");
    const scalar energyJump = offsets.get<scalar>(vapourName)-offsets.get<scalar>(liquidName);
    const scalar dt = mesh_.time().deltaTValue();
    const volScalarField rhoLiquid(liquid->thermo().rho());
    const volScalarField rhoVapour(vapour->thermo().rho());
    const volScalarField cvLiquid(liquid->thermo().Cv());
    const volScalarField cvVapour(vapour->thermo().Cv());
    const scalar cvMin = min(gMin(cvLiquid.primitiveField()), gMin(cvVapour.primitiveField()));
    const scalar cvMax = max(gMax(cvLiquid.primitiveField()), gMax(cvVapour.primitiveField()));
    if (!std::isfinite(cvMin) || !std::isfinite(cvMax)
     || cvMin <= 0 || cvMax-cvMin > 1e-10*cvMax)
    {
        FatalErrorInFunction << "Interface energy candidate requires equal constant Cv"
            << exit(FatalError);
    }
    const volScalarField molWeight(liquid->thermo().W());
    const volScalarField areaDensity(mag(fvc::grad(*liquid)));
    scalar boilingTemperature = 0, latentHeat = 0;
    const auto boiling = boils_.find(interfacePair(*liquid, *vapour));
    const auto latent = LatentHeatGass_.find(interfacePair(*liquid, *vapour));
    if (boiling == boils_.end() || latent == LatentHeatGass_.end())
    {
        FatalErrorInFunction << "Missing boiling/latent data for interface pair" << exit(FatalError);
    }
    boilingTemperature = boiling();
    latentHeat = latent();
    const dimensionedScalar referencePressure("P0", dimPressure, phasedictionary);
    scalar area = 0, requestedRate = 0, transferredRate = 0, limitedRate = 0;
    scalar pairResidual = 0, volumeRate = 0, power = 0;
    forAll(areaDensity, celli)
    {
        const scalar rl = rhoLiquid[celli], rv = rhoVapour[celli], temperature = T_[celli];
        if (rl <= 0 || rv <= 0 || temperature <= 0 || dt <= 0 || boilingTemperature <= 0)
        {
            FatalErrorInFunction << "Non-positive interface thermodynamic state" << exit(FatalError);
        }
        const scalar a = max(scalar(0), min(scalar(1), (*liquid)[celli]));
        const scalar b = max(scalar(0), min(scalar(1), (*vapour)[celli]));
        scalar j = prescribedJ;
        if (interfaceModel == "kinetic")
        {
            // Pure vapour only: p is backpressure; not a shielding-gas partial pressure model.
            const scalar molarMass = molWeight[celli]/1000;
            const scalar psat = referencePressure.value()*Foam::exp
            (
                molarMass*latentHeat/(8.314*boilingTemperature)
               *(1-boilingTemperature/temperature)
            );
            j = accommodation*(psat-p_[celli])
               *Foam::sqrt(molarMass/(2*M_PI*8.314*temperature));
        }
        if (!std::isfinite(j) || !std::isfinite(energyJump) || !std::isfinite(rl)
         || !std::isfinite(rv) || !std::isfinite(temperature))
        {
            FatalErrorInFunction << "Non-finite interface state/flux" << exit(FatalError);
        }
        const scalar requested = scalar(phaseChangeEnabled)*j*areaDensity[celli];
        // Bounds include the volume-expansion correction used by explicit alpha.
        const scalar cLiquid = (1-a)/rl+a/rv;
        const scalar cVapour = (1-b)/rv+b/rl;
        const scalar positiveLimit = inventoryFraction/dt
           *min(a/cLiquid, (1-b)/cVapour);
        const scalar negativeLimit = inventoryFraction/dt
           *min(b/cVapour, (1-a)/cLiquid);
        const scalar transferred = min(positiveLimit, max(-negativeLimit, requested));
        interfaceAlphaSources[liquidIndex][celli] = -transferred/rl;
        interfaceAlphaSources[vapourIndex][celli] = transferred/rv;
        commonPhaseChangePower[celli] = energyJump*transferred;
        const scalar v = mesh_.V()[celli];
        area += areaDensity[celli]*v;
        requestedRate += requested*v;
        transferredRate += transferred*v;
        limitedRate += mag(requested-transferred)*v;
        pairResidual += (rl*interfaceAlphaSources[liquidIndex][celli]
                       +rv*interfaceAlphaSources[vapourIndex][celli])*v;
        volumeRate += transferred*(1/rv-1/rl)*v;
        power += commonPhaseChangePower[celli]*v;
    }
    for (scalar* value : {&area, &requestedRate, &transferredRate, &limitedRate, &pairResidual, &volumeRate, &power})
    {
        reduce(*value, sumOp<scalar>());
    }
    if (Pstream::master())
    {
        const int savedPrecision = Info().precision();
        Info().precision(15);
        Info<< "M1C_SOURCE time=" << mesh_.time().value() << " dt=" << dt
            << " model=" << interfaceModel << " enabled=" << phaseChangeEnabled
            << " areaM2=" << area << " requestedKgPerS=" << requestedRate
            << " transferredKgPerS=" << transferredRate << " limitedAbsKgPerS=" << limitedRate
            << " pairMassResidualKgPerS=" << pairResidual << " volumeSourceM3PerS=" << volumeRate
            << " phasePowerW=" << power << " energyJumpJPerKg=" << energyJump << nl;
        Info().precision(savedPrecision);
    }
}
