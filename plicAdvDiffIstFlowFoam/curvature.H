
alpha1_smooth = fvc::average(alpha1);
alpha1_smooth = fvc::average(alpha1_smooth);
alpha1_smooth = fvc::average(alpha1_smooth);

const surfaceVectorField& Sf = mesh.Sf();

// Cell gradient of alpha
const volVectorField gradAlpha(fvc::grad(alpha1_smooth, "nHat"));

// Interpolated face-gradient of alpha
surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

// Face unit interface normal
surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN));

// Face unit interface normal flux
nHatf = nHatfv & Sf;

// Simple expression for curvature
K = -fvc::div(nHatf);
