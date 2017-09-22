
if (adjustTimeStep)
{
    if ((runTime.timeIndex() == 0) && (diffCoNum > SMALL))
    {
        runTime.setDeltaT
        (
            min
            (
                maxDiffCo*runTime.deltaTValue()/diffCoNum,
                runTime.deltaTValue()
            )
        );
    }
}
