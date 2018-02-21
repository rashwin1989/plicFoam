
if (adjustTimeStep)
{
    if ((runTime.timeIndex() == 0))
    {
        scalar maxDeltaTFactAdv;
        scalar maxDeltaTFactIFT;
        scalar maxDeltaTFact;

        if(CoNum > SMALL)
        {
            maxDeltaTFactAdv = maxCo/CoNum;
        }
        else
        {
            maxDeltaTFactAdv = 1;
        }        

        if(IFTCoNum > SMALL)
        {
            maxDeltaTFactIFT = maxIFTCo/IFTCoNum;
        }
        else
        {
            maxDeltaTFactIFT = 1;
        }

        maxDeltaTFact = min(maxDeltaTFactAdv, maxDeltaTFactIFT);

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
