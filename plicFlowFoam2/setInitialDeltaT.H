
if (adjustTimeStep)
{
    if ((runTime.timeIndex() == 0))
    {
        scalar maxDeltaTFactAdv;        
        scalar maxDeltaTFact;

        if(CoNum > SMALL)
        {
            maxDeltaTFactAdv = maxCo/CoNum;
        }
        else
        {
            maxDeltaTFactAdv = 1;
        }        

        maxDeltaTFact = maxDeltaTFactAdv;

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
