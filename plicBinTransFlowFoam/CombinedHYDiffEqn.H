
if(diff_debug)
{
    plicFuncs::print_line(osDiff, 80);
    osDiff<< "Entered CombinedHYDiffEqn.H" << nl
        << "Solution Time = " << runTime.timeName() << "  Y Diff" << endl;
    plicFuncs::print_line(osDiff, 80);
    osDiff<< endl;
}

//set fractional time-step for diffusion
deltaT_diff = deltaT;
dt = deltaT_diff;

#include "correct_thermo_trans_prop.H"
#include "diff_grad_interp.H" 

if(diff_debug)
{
    plicFuncs::print_line(osDiff, 80);
    osDiff<< "Done first gradienct calculation" << nl
        << "Solution Time = " << runTime.timeName() << "  Y Diff" << endl;
    plicFuncs::print_line(osDiff, 80);
    osDiff<< endl;
}

#include "HCondDiffEqn.H"

if(diff_debug)
{
    plicFuncs::print_line(osDiff, 80);
    osDiff<< "Done first half HCondDiff.H" << nl
        << "Solution Time = " << runTime.timeName() << "  Y Diff" << endl;
    plicFuncs::print_line(osDiff, 80);
    osDiff<< endl;
}

#include "YHDiffEqn.H"

if(diff_debug)
{
    plicFuncs::print_line(osDiff, 80);
    osDiff<< "Done YHDiffEqn.H" << nl
        << "Solution Time = " << runTime.timeName() << "  Y Diff" << endl;
    plicFuncs::print_line(osDiff, 80);
    osDiff<< endl;
}

#include "correct_thermo_trans_prop.H"
#include "diff_grad_interp.H" 

if(diff_debug)
{
    plicFuncs::print_line(osDiff, 80);
    osDiff<< "Done second  gradienct calculation" << nl
        << "Solution Time = " << runTime.timeName() << "  Y Diff" << endl;
    plicFuncs::print_line(osDiff, 80);
    osDiff<< endl;
}

#include "HCondDiffEqn.H"

if(diff_debug)
{
    plicFuncs::print_line(osDiff, 80);
    osDiff<< "Done second half HCondDiff.H" << nl
        << "Solution Time = " << runTime.timeName() << "  Y Diff" << endl;
    plicFuncs::print_line(osDiff, 80);
    osDiff<< endl;
}

#include "correct_thermo_trans_prop.H"
#include "diff_grad_interp.H" // For UEqn.H that used mu1f
