
if (adjustTimeStep)
{
    if ((runTime.timeIndex() == 0))
    {
        scalar maxDeltaTFactAdv;
        scalar maxDeltaTFactDiff;
        scalar maxDeltaTFact;

        if(CoNum > SMALL)
        {
            maxDeltaTFactAdv = maxCo/CoNum;
        }
        else
        {
            maxDeltaTFactAdv = 1;
        }

        if(diffCoNum > SMALL)
        {
            maxDeltaTFactDiff = maxDiffCo/diffCoNum;
        }
        else
        {
            maxDeltaTFactDiff = 1;
        }

        maxDeltaTFact = min(min(maxDeltaTFactAdv, maxDeltaTFactDiff), 5);

        runTime.setDeltaT
        (
            min
            (
                maxDeltaTFact*runTime.deltaTValue(),
                maxDeltaT
            )
        );

        Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
    }
}
