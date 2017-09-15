
scalar t = runTime.value();

Info<< "Setting shear velocity field for time = " << t << endl;

vectorField& Ui = U.internalField();
const vectorField& Ci = mesh.C().internalField();

scalar pi = Foam::constant::mathematical::pi;

//bool twoD = false;

forAll(Ui, cellI)
{
    if(twoD)
    {
        scalar x = Ci[cellI].x() + 0.5;
        scalar y = Ci[cellI].y() + 0.5;

        scalar sinx = Foam::sin(pi*x);
        scalar cosx = Foam::cos(pi*x);
        scalar siny = Foam::sin(pi*y);
        scalar cosy = Foam::cos(pi*y);
        scalar cost = Foam::cos(pi*t/8);

        scalar Ux = -2*sinx*sinx*siny*cosy*cost;
        scalar Uy = 2*siny*siny*sinx*cosx*cost;

        Ui[cellI].x() = Ux;
        Ui[cellI].y() = Uy;
    }
    else
    {
        scalar x = Ci[cellI].x();
        scalar y = Ci[cellI].y();        
        scalar z = Ci[cellI].z();

        scalar sinx = Foam::sin(pi*x);
        scalar siny = Foam::sin(pi*y);
        scalar sinz = Foam::sin(pi*z);
        scalar sin2x = Foam::sin(2*pi*x);
        scalar sin2y = Foam::sin(2*pi*y);
        scalar sin2z = Foam::sin(2*pi*z);
        scalar cost = Foam::cos(pi*t/3);

        scalar Ux = 2*sinx*sinx*sin2y*sin2z*cost;
        scalar Uy = -sin2x*siny*siny*sin2z*cost;
        scalar Uz = -sin2x*sin2y*sinz*sinz*cost;

        Ui[cellI].x() = Ux;
        Ui[cellI].y() = Uy;
        Ui[cellI].z() = Uz;
    }
}

U.correctBoundaryConditions();

phi = linearInterpolate(U) & mesh.Sf();
