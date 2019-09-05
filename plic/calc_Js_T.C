void calc_Js_T
(
    const fvMesh& mesh,
    const labelListList& cellStencil,    
    const List<List<scalar> >& Y1_flatFld,
    const List<List<scalar> >& Y0_flatFld,
    const List<scalar>& T1_flatFld,
    const List<scalar>& T0_flatFld,
    const List<scalar>& alpha1_flatFld,
    const volScalarField& alpha1,
    const volScalarField& rho1,
    const volScalarField& rho0,
    const PtrList<volScalarField>& D1,
    const PtrList<volScalarField>& D0,
    const PtrList<volScalarField>& hpar1,
    const PtrList<volScalarField>& hpar0,
    const List<vector>& C_ph1_flatFld,
    const List<vector>& C_ph0_flatFld,
    const volVectorField& C_intfc,
    const volScalarField& A_intfc,
    const volVectorField& nHat,
    const label& nSpecies,
    const volScalarField& T1, 
    const volScalarField& T0, 
    const volScalarField& Ts, 
    const scalar& ALPHA_2PH_MIN,
    double erf_a,
    double erf_b,
    int nErf,
    double *erfInv_table,
    bool useErf,
    const List<scalar>& Yinf1,
    const List<scalar>& Yinf0,
    const PtrList<volScalarField>& Ys1,
    const PtrList<volScalarField>& Ys0,
    PtrList<volScalarField>& Hs1,
    PtrList<volScalarField>& Hs0,
    PtrList<volScalarField>& grads1,
    PtrList<volScalarField>& grads0,
    volScalarField& delta1,
    volScalarField& delta0,
    const bool debug,
    OFstream& os
)
{
    int n = nSpecies;
    int i;
    double erf1, x_delta1, delta1_cellI, erf0, x_delta0, delta0_cellI, pi, two_sqrtPi;
    scalar alpha1_cellI, A_intfc_cellI, dn1, dn0, Teff1, Teff0;
    vector nf, C_intfc_cellI;

    List<scalar> grads1_cellI(n);
    List<scalar> grads0_cellI(n);
    List<scalar> Yeff1(n);
    List<scalar> Yeff0(n);
    List<scalar> ys1(n);
    List<scalar> ys0(n);

    pi = Foam::constant::mathematical::pi;
    two_sqrtPi = 2.0/sqrt(pi);

    scalar ALPHA_2PH_MAX = 1 - ALPHA_2PH_MIN;

    //const labelList& own = mesh.owner();
    //const labelList& nei = mesh.neighbour();

    //const surfaceVectorField& Sf = mesh.Sf();
    //const surfaceScalarField& magSf = mesh.magSf();
    //const surfaceVectorField& Cf = mesh.Cf();

    const scalarField& alpha1Cells = alpha1.internalField();    
    const vectorField& C_intfcCells = C_intfc.internalField();
    const scalarField& A_intfcCells = A_intfc.internalField();
    const vectorField& nHatCells = nHat.internalField();
    const vectorField& T1Cells = T1.internalField();
    const vectorField& T0Cells = T0.internalField();
    scalarField& TsCells = Ts.internalField();
    scalarField& delta1Cells = delta1.internalField();
    scalarField& delta0Cells = delta0.internalField();

    if(debug)
    {
        print_line(os, 100);
        os<< "Interfacial Heat Conduction Calculation" << endl;
        print_line(os, 100);
        os<< endl;
        print_line(os, 100);
        os<< "Internal cells" << endl;
        print_line(os, 100);
        os<< endl;
    }

    //Hs for all interface cells
    forAll(alpha1Cells, cellI)
    {
        alpha1_cellI = alpha1Cells[cellI];        

        double Ts_cellI = alpha1_cellI*T1Cells[cellI] + (1 - alpha1_cellI)*T0Cells[cellI];

        if(alpha1_cellI > ALPHA_2PH_MIN && alpha1_cellI < ALPHA_2PH_MAX)
        {
            nf = nHatCells[cellI];
            C_intfc_cellI = C_intfcCells[cellI];
            A_intfc_cellI = A_intfcCells[cellI];
            labelList curCellsAll = cellStencil[cellI];

            /*            for(i=0; i<n; i++)
            {
                ys1[i] = Ys1[i].internalField()[cellI];
                ys0[i] = Ys0[i].internalField()[cellI];
            } */

            if(debug)
            {
                print_line(os, 100);
                os<< "Cell: " << cellI << "  alpha1 = " << alpha1_cellI << endl;
                print_line(os, 100);
                os<< "nf = " << nf << "  C_intfc = " << C_intfc_cellI << "  A_intfc = " << A_intfc_cellI << nl
                    << "C_ph1 = " << C_ph1_flatFld[cellI] << "  C_ph0 = " << C_ph0_flatFld[cellI] << nl
                    << "rho1 = " << rho1Cells[cellI] << "  rho0 = " << rho0Cells[cellI] << endl;
                print_line(os, 100);
                os<< endl;
            }            

            //phase-1
            //ensure nf direction is into the phase
            //            calc_cell_intfcGrad_coeffs(mesh, cellI, nf, C_intfc_cellI, Y1_flatFld, alpha1_flatFld, C_ph1_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 1, dn1, Yeff1, debug, os);
            calc_cell_intfcGrad_coeffs(mesh, cellI, nf, C_intfc_cellI, Y1_flatFld, T1_flatFld, alpha1_flatFld, C_ph1_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 1, dn1, Yeff1, Teff1, debug, os);

            if(useErf)
            {            
            //Calculate the BL thickness delta from error function profile for Y[0]            
            erf1 = (Yeff1[0] - ys1[0])/(Yinf1[0] - ys1[0]);
            erf1 = max(erf1, 0.0);
            erf1 = min(erf1, 1.0);
            calc_erfInv_from_table(erf1,erf_a,erf_b,nErf,erfInv_table,x_delta1);
            delta1_cellI = dn1/x_delta1;

            delta1Cells[cellI] = delta1_cellI;
            }
            
            //phase-0
            //ensure nf direction is into the phase
            //then reverse nf again for Js0 calculation            
            //            calc_cell_intfcGrad_coeffs(mesh, cellI, -nf, C_intfc_cellI, Y0_flatFld, alpha1_flatFld, C_ph0_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 0, dn0, Yeff0, debug, os);
            calc_cell_intfcGrad_coeffs(mesh, cellI, -nf, C_intfc_cellI, Y0_flatFld, T0_flatFld, alpha1_flatFld, C_ph0_flatFld, curCellsAll, nSpecies, ALPHA_2PH_MIN, 0, dn0, Yeff0, Teff0, debug, os);

            if(useErf)
            {            
            //Calculate the BL thickness delta from error function profile for Y[0]            
            erf0 = (Yeff0[0] - ys0[0])/(Yinf0[0] - ys0[0]);
            erf0 = max(erf0, 0.0);
            erf0 = min(erf0, 1.0);
            calc_erfInv_from_table(erf0,erf_a,erf_b,nErf,erfInv_table,x_delta0);
            delta0_cellI = dn0/x_delta0;

            delta0Cells[cellI] = delta0_cellI;
            }

            scalar lambda1_cellI = lambda1.internalField()[cellI];
            scalar lambda0_cellI = lambda0.internalField()[cellI];
            Hs1.internalField[cellI] =  A_intfc_cellI*lambda1_cellI*(Ts_cellI - Teff1)/dn1;
            Hs0.internalField[cellI] = - Hs1.internalField[cellI];
            TsCells[cellI] = (lambda1_cellI*Teff1/dn1 + lambda0_cellI*Teff0/dn0)/(lambda1_cellI/dn1 + lambda0_cellI/dn0);

            if(debug && useErf)
            {
                print_line(os, 100);
                os<< "Phase-1:" << endl;
                os<< "erf1 = " << erf1 << "  x_delta1 = " << x_delta1 << "  dn1 = " << dn1 << "  delta1 = " << delta1_cellI << endl;
                print_line(os, 100);
                os<< setw(9) << "Species" << setw(14) << "ys1" << setw(14) << "Yeff1" << setw(14) << "Yinf1" << setw(14) << "grads1" << setw(14) << "D1" << setw(14) << "Js1" << endl;
                print_line(os, 100);
                for(i=0; i<n; i++)
                {
                    os<< setw(9) << i << setw(14) << ys1[i] << setw(14) << Yeff1[i] << setw(14) << Yinf1[i] << setw(14) << grads1_cellI[i] << setw(14) << D1[i].internalField()[cellI] << setw(14) << Js1[i].internalField()[cellI] << endl;
                }
                print_line(os, 100);

                print_line(os, 100);
                os<< "Phase-0:" << endl;
                os<< "erf0 = " << erf0 << "  x_delta0 = " << x_delta0 << "  dn0 = " << dn0 << "  delta0 = " << delta0_cellI << endl;
                print_line(os, 100);
                os<< setw(9) << "Species" << setw(14) << "ys0" << setw(14) << "Yeff0" << setw(14) << "Yinf0" << setw(14) << "grads0" << setw(14) << "D0" << setw(14) << "Js0" << endl;
                print_line(os, 100);
                for(i=0; i<n; i++)
                {
                    os<< setw(9) << i << setw(14) << ys0[i] << setw(14) << Yeff0[i] << setw(14) << Yinf0[i] << setw(14) << grads0_cellI[i] << setw(14) << D0[i].internalField()[cellI] << setw(14) << Js0[i].internalField()[cellI] << endl;
                }
                print_line(os, 100);
                os<< endl;
            }
        }
        else
        {
            for(i=0; i<nSpecies; i++)
            { 
                Hs1.internalField()[cellI] = 0;
                Hs0.internalField()[cellI] = 0;
                TsCells[cellI] = Ts_cellI;
            }
        }        
    }

    if(debug)
    {
        print_line(os, 100);
        os<< " Done Interfacial Conduction Calculation" << endl;
        print_line(os, 100);
        os<< endl;
    }
}

