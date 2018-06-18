
plicFuncs::print_line(osU, 100);
osU<< "Solution Time = " << runTime.timeName() << "  Momentum Eqn" << endl;
plicFuncs::print_line(osU, 100);
osU<< endl;

alpha1f = fvc::interpolate(limitedAlpha1, phi, "alpha1f");
alpha0f = scalar(1) - alpha1f;
muf = alpha1f*mu1f + alpha0f*mu0f;

if(U_debug)
{
    plicFuncs::print_field(mu1, osU); 
    plicFuncs::print_field(mu0, osU);
    plicFuncs::print_surfaceField(alpha1f, mesh, osU);
    plicFuncs::print_surfaceField(mu1f, mesh, osU);
    plicFuncs::print_surfaceField(mu0f, mesh, osU);
    plicFuncs::print_surfaceField(muf, mesh, osU);
}

/*
surfaceScalarField gammaf = mag(interface.nHatf())/mesh.magSf();
gammaf.min(1.0);
surfaceScalarField muf_a = alpha1f*mu1f + alpha0f*mu0f;
surfaceScalarField muf_h = mu1f*mu0f/muf_a;
muf = gammaf*muf_h + (1 - gammaf)*muf_a;
*/

Info<< "Min(rho) = " << min(rho).value()
    << "  Max(rho) = " << max(rho).value()
    << endl;
Info<< "Min(mu) = " << min(mu).value()
    << "  Max(mu) = " << max(mu).value()
    << endl;
Info<< "Min(muf) = " << min(muf).value()
    << "  Max(muf) = " << max(muf).value()
    << nl
    << endl;


fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
    + fvm::div(rhoPhi, U)
    - fvm::laplacian(muf, U)
);

UEqn.relax();

if (momentumPredictor)
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
                fvc::interpolate(sigma*K)*fvc::snGrad(alpha1)
              - ghf*fvc::snGrad(rho)
              - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        )
    );
}
