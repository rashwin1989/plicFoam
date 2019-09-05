void calc_mS_Hs_alphaS
(
    const fvMesh& mesh,
    const PtrList<volScalarField>& C1,
    const PtrList<volScalarField>& C0,
    const PtrList<volScalarField>& Y1,
    const PtrList<volScalarField>& Y0,
    const volScalarField& T1,
    const volScalarField& T0,
    const volScalarField& alpha1,
    const volScalarField& rho1,
    const volScalarField& rho0,
    const PtrList<volScalarField>& Ys1,
    const PtrList<volScalarField>& Ys0,
    const PtrList<volScalarField>& Js1,
    const PtrList<volScalarField>& Js0,
    const PtrList<volScalarField>& hpar1,
    const PtrList<volScalarField>& hpar0,
    const label& nSpecies,
    const scalar& ALPHA_2PH_MIN,
    PtrList<volScalarField>& mS1,
    PtrList<volScalarField>& mS0,
    volScalarField& mS1Tot,
    volScalarField& mS0Tot,
    volScalarField& Hs1,
    volScalarField& Hs0,
    volScalarField& alphaS1,
    volScalarField& alphaS0,
    const scalar& dt,
    const bool debug,
    OFstream& os
)
{
    int n, i;
    n = nSpecies;
    scalar alpha1_cellI, rho1_cellI, rho0_cellI, V_cellI, mS1Tot_cellI, mS1Tot_cellI_tmp, limiterTot, limiter_min;

    List<scalar> mS1_cellI(n);
    List<scalar> Js1_cellI(n);
    List<scalar> Js0_cellI(n);
    List<scalar> Ys1_cellI(n);
    List<scalar> Ys0_cellI(n);
    List<scalar> Y1_cellI(n);
    List<scalar> Y0_cellI(n);
    List<scalar> C1_cellI(n);
    List<scalar> C0_cellI(n);
    List<scalar> limiterY(n);

    scalar Hs1_cellI;    
    scalar ALPHA_2PH_MAX = 1 - ALPHA_2PH_MIN;
    const scalarField& V = mesh.V();

    const scalarField& alpha1Cells = alpha1.internalField();    
    const scalarField& rho1Cells = rho1.internalField();
    const scalarField& rho0Cells = rho0.internalField();
    const scalarField& hpar1Cells = hpar1.internalField();
    const scalarField& hpar0Cells = hpar0.internalField();
    scalarField& mS1TotCells = mS1Tot.internalField();
    scalarField& mS0TotCells = mS0Tot.internalField();
    scalarField& alphaS1Cells = alphaS1.internalField();
    scalarField& alphaS0Cells = alphaS0.internalField();
    scalarField& Hs1cells = Hs1.internalField();
    scalarField& Hs0cells = Hs0.internalField();

    if(debug)
    {
        print_line(os, 100);
        print_line(os, 100);
        os<< "Interfacial Species Transfer Source Terms Calculation" << endl;
        print_line(os, 100);
        print_line(os, 100);
        os<< "Internal cells" << endl;
        print_line(os, 100);
    }

    //mS for all interface cells
    forAll(alpha1Cells, cellI)
    {
        alpha1_cellI = alpha1Cells[cellI];        

        if(alpha1_cellI > ALPHA_2PH_MIN && alpha1_cellI < ALPHA_2PH_MAX)
        {
            rho1_cellI = rho1Cells[cellI];        
            rho0_cellI = rho0Cells[cellI];        
            V_cellI = V[cellI];        
            for(i=0; i<n; i++)
            {
                Js1_cellI[i] = Js1[i].internalField()[cellI]/V_cellI;
                Js0_cellI[i] = Js0[i].internalField()[cellI]/V_cellI;
                Ys1_cellI[i] = Ys1[i].internalField()[cellI];
                Ys0_cellI[i] = Ys0[i].internalField()[cellI];
                Y1_cellI[i] = Y1[i].internalField()[cellI];
                Y0_cellI[i] = Y0[i].internalField()[cellI];
                C1_cellI[i] = C1[i].internalField()[cellI];
                C0_cellI[i] = C0[i].internalField()[cellI];
            }
            mS1Tot_cellI = 0; mS1Tot_cellI_tmp = 0;
            
            mS1Tot_cellI_tmp = Js0_cellI[0] - Js1_cellI[0];
            if(mS1Tot_cellI_tmp > 0)
            {
                if((Ys1_cellI[0] - Y0_cellI[0]) > 0)
                {
                    mS1Tot_cellI_tmp /= (Ys1_cellI[0] - Y0_cellI[0]);
                }
                else
                {
                    mS1Tot_cellI_tmp /= (Ys1_cellI[0] - Ys0_cellI[0]);
                }

                for(i=0; i<n; i++)
                {                                
                    mS1_cellI[i] = mS1Tot_cellI_tmp*Ys1_cellI[i] + Js1_cellI[i];
                }
            }
            else
            {
                if((Y1_cellI[0] - Ys0_cellI[0]) > 0)
                {
                    mS1Tot_cellI_tmp /= (Y1_cellI[0] - Ys0_cellI[0]);
                }
                else
                {
                    mS1Tot_cellI_tmp /= (Ys1_cellI[0] - Ys0_cellI[0]);
                }

                for(i=0; i<n; i++)
                {                                
                    mS1_cellI[i] = mS1Tot_cellI_tmp*Ys0_cellI[i] + Js0_cellI[i];
                }
            }

            if(debug)
            {
                print_line(os, 100);
                print_line(os, 100);
                os<< "Cell: " << cellI << "  alpha1 = " << alpha1_cellI << endl;
                print_line(os, 100);
                os<< "mS limiter calculation" << endl;
                print_line(os, 100);
            }

            limiterTot = 1;
            limiter_min = 1;
            for(i=0; i<n; i++) limiterY[i] = 1;

            calc_mS_limiter(C1_cellI, C0_cellI, Y1_cellI, Y0_cellI, mS1_cellI, mS1Tot_cellI_tmp, dt, n, limiterTot, limiterY, limiter_min, debug, os);     
     
            for(i=0; i<n; i++)
            {                                
                mS1_cellI[i] = limiter_min*mS1_cellI[i];                

                mS1[i].internalField()[cellI] = mS1_cellI[i];
                mS0[i].internalField()[cellI] = -mS1_cellI[i];
            }
            
            mS1Tot_cellI = limiter_min*mS1Tot_cellI_tmp;

            mS1TotCells[cellI] = mS1Tot_cellI;
            mS0TotCells[cellI] = -mS1Tot_cellI;

            alphaS1Cells[cellI] = mS1Tot_cellI/rho1_cellI;
            alphaS0Cells[cellI] = -mS1Tot_cellI/rho0_cellI;

            Hs1_cellI = 0;
            for(i=0; i<n; i++)
            {
                scalar mSi = mS1[i].internalField()[cellI];
                Hs1_cellI += max(mSi,0)*hpar0.internalField()[cellI] + min(mSi,0)*hpar1.internalField()[cellI];
            }
            Hs1.internalField()[cellI] = Hs1_cellI;
            Hs0.internalField()[cellI] = -Hs1_cellI;

            if(debug)
            {
                print_line(os, 100);
                print_line(os, 100);
                os<< "Cell: " << cellI << "  alpha1 = " << alpha1_cellI << endl;
                print_line(os, 100);
                os<< setw(8) << "Species" << "  " << setw(14) << "C1" << "  " << setw(14) << "C0" << "  " << setw(14) << "Y1" << "  " << setw(14) << "Y0" << endl;
                print_line(os, 100);
                for(i=0; i<n; i++)
                {
                    os<< setw(7) << i << "  " << setw(14) << C1[i].internalField()[cellI] << "  " << setw(14) << C0[i].internalField()[cellI] << "  " << setw(14) << Y1[i].internalField()[cellI] << "  " << setw(14) << Y0[i].internalField()[cellI] << endl;
                }
                print_line(os, 100);
                os<< setw(8) << "Species" << "  " << setw(14) << "Js1" << "  " << setw(14) << "Js0" << "  " << setw(14) << "mS1" << "  " << setw(14) << "mS0" << endl;
                print_line(os, 100);
                for(i=0; i<n; i++)
                {
                    os<< setw(7) << i << "  " << setw(14) << Js1[i].internalField()[cellI] << "  " << setw(14) << Js0[i].internalField()[cellI] << "  " << setw(14) << mS1[i].internalField()[cellI] << "  " << setw(14) << mS0[i].internalField()[cellI] << endl;
                }
                print_line(os, 100);
                os<< "limiter_min = " << limiter_min << nl
                    << "mS1Tot = " << mS1TotCells[cellI] << "  mS0Tot = " << mS0TotCells[cellI] << nl
                    << "rho1 = " << rho1_cellI << "  rho0 = " << rho0_cellI << nl
                    << "alphaS1 = " << alphaS1Cells[cellI] << "  alphaS0 = " << alphaS0Cells[cellI] << endl;
                print_line(os, 100);
                print_line(os, 100);
            }
        }
        else
        {
            for(i=0; i<n; i++)
            {
                mS1[i].internalField()[cellI] = 0;
                mS0[i].internalField()[cellI] = 0;
            }

            mS1TotCells[cellI] = 0;
            mS0TotCells[cellI] = 0;

            alphaS1Cells[cellI] = 0;
            alphaS0Cells[cellI] = 0;

            Hs1Cells[cellI] = 0;
            Hs0Cells[cellI] = 0;
        }        
    }

    if(debug)
    {
        print_line(os, 100);
        print_line(os, 100);
        os<< "Done Interfacial Species Transfer Source Terms Calculation" << endl;
        print_line(os, 100);
        print_line(os, 100);
    }
}
