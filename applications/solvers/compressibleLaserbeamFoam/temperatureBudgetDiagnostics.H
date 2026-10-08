// Observe the unrelaxed rho*T equation at the final thermal correction.
// These are integrated temperature-equation terms (kg K/s), not a general
// thermodynamic energy balance. Multiplication by Cv is valid only for the
// uniform, common constant-Cv synthetic fixture.
if (runTime.controlDict().lookupOrDefault<Switch>("temperatureBudgetDiagnostics", false))
{
    const scalar storage = gSum
    (
        fvc::ddt(rho, T)().primitiveField()*mesh.V().field()
    );
    const scalar transport = gSum
    (
        fvc::div(mixture.rhoPhi(), T)().primitiveField()*mesh.V().field()
    );
    const scalar diffusion = -gSum
    (
        fvc::laplacian(mixture.alphaEff(turbulence->mut()), T)().primitiveField()
       *mesh.V().field()
    );
    const tmp<volScalarField> mechanicalTerm =
    (
        fvc::div(fvc::absolute(phi, U), p)
      + fvc::ddt(rho, K) + fvc::div(mixture.rhoPhi(), K)
    )*mixture.rCv();
    const scalar mechanical = gSum
    (
        mechanicalTerm().primitiveField()*mesh.V().field()
    );
    const tmp<volScalarField> laserTerm = mixture.rCv()*laser.deposition();
    const tmp<volScalarField> fusionTerm = mixture.rCv()*TRHS;
    const scalar laserSource = -gSum(laserTerm().primitiveField()*mesh.V().field());
    const scalar fusion = gSum(fusionTerm().primitiveField()*mesh.V().field());
    const scalar phaseChange = gSum(mass_dot.primitiveField()*mesh.V().field());
    const scalar inventory = gSum
    (
        rho.primitiveField()*T.primitiveField()*mesh.V().field()
    );
    const scalar oldInventory = gSum
    (
        rho.oldTime().primitiveField()*T.oldTime().primitiveField()*mesh.V().field()
    );
    if (Pstream::master())
    {
        const int oldPrecision = Info().precision();
        Info().precision(15);
        Info<< "M1_T_BUDGET time=" << runTime.value()
            << " thermalCorrections=" << iter
            << " inventoryKgK=" << inventory
            << " oldInventoryKgK=" << oldInventory
            << " storageKgKPerS=" << storage
            << " transportKgKPerS=" << transport
            << " diffusionKgKPerS=" << diffusion
            << " mechanicalKgKPerS=" << mechanical
            << " laserKgKPerS=" << laserSource
            << " fusionKgKPerS=" << fusion
            << " phaseChangeKgKPerS=" << phaseChange
            << " residualKgKPerS="
            << storage+transport+diffusion+mechanical+laserSource+fusion+phaseChange
            << nl;
        Info().precision(oldPrecision);
    }
}
