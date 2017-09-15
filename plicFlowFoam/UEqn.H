
surfaceScalarField alpha1f = fvc::interpolate(limitedAlpha1, phi, "alpha1f");
surfaceScalarField alpha0f = scalar(1) - alpha1f;
surfaceScalarField gammaf = mag(interface.nHatf())/mesh.magSf();
gammaf.min(1.0);
surfaceScalarField muf_a = alpha1f*mu1 + alpha0f*mu0;
surfaceScalarField muf_h = mu1*mu0/muf_a;
muf = gammaf*muf_h + (1 - gammaf)*muf_a;

Info<< "Min(rho) = " << min(rho).value()
    << "  Max(rho) = " << max(rho).value()
    << endl;
Info<< "Min(mu) = " << min(mu).value()
    << "  Max(mu) = " << max(mu).value()
    << endl;
Info<< "Min(muf) = " << min(muf).value()
    << "  Max(muf) = " << max(muf).value()
    << endl;

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
    + fvm::div(rhoPhi, U)
    - fvm::laplacian(muf, U)
);

UEqn.relax();

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
                fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
              - ghf*fvc::snGrad(rho)
              - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        )
    );
}
