# Graph Learning and Augmentation Based Interpolation of Signal Strength for Location-aware Communicationsr

## Intorduction
The MATLAB programs contain in this directory solving the power map interpolation problem using both Graph Learning and Augmentation [1], Gaussian Process [2] and Kriging [3], as well as functions that perform Monta Carlo simulation.

## CONTANTS
In main directory :
1.GLA_interpolator.m : main function for GLA interpolator [1]
2.GP_interpolator.m : main function for GP interpolator [2]
3.LS_interpolator.m : main function for Kriging interpolator [3]
4.AggregatedMSETest.m : Monta Carlo simulation for aggregated MSE
5.RegionalMSETest.m : Monta Carlo simulation for Regonal MSE
6.SNRTest.m : Monta Carlo simulation for SNR v.s.MSE
7.startup.m : add path and make directory, executed automatically at startup
8.cleanup.m : clean all *.mat file in Data/ directory
    In Function/ directory : 
        MakeGroundTrue.m : generate ground true power map [4]
        MakeObservation.m : generate observed power map
        MakeSample.m : generate sample points
        DataStatistic.m : generate data statistic
        GMRF.m : GMRF graph learning [5]
        LearnParameter.m : learn paramters for graph augmentation
    In BCD_Algorithm/ diretory [5] : 
        EstimateGGL.m : GGL graph learning
        EstimateDDGL.m : DDGL graph learning
        EstimateCGL.m : CGL graph learning
        nonnegative_qp_solver.m : function for solving inner QP subproblem
        update_sherman_morrison_diag.m : funciton for updating diagonal elements in DDGL graph learning
    In MonteCarlo/ directory : 
        AggregatedMSETest/ : store outputs of AggregatedMSETest.m
        RegionalMSETest/ :store outputs of RegionalMSETest.m
        SNRTest/ :store outputs of SNRTest.m
        PlotAggregatedMSETest.m : plot results of AggregatedMSETest.m
        PlotRegionalMSETest.m : plot results of RegionalMSETest.m
        PlotSNRTest.m : plot results of SNRTest.m
    In Data/ directory :
        GP/ : store results of GP_interpolator.m
        GraphLearningData/ : store results of DataStatistic.m, GMRF.m and LearnParameter.m
        GroundTrueData/ : store results of MakeGroundTrue.m
        Sample/ : store results of MakeSample.m
    In Figure/ directory :
        TruePowerMap.fig : true power map generated from MakeGroundTrue.m
        ObservedPowerMap.fig : noisy observation of power map plotted in TruePowerMap.fig
        InterpolatedPowerMap(MAP).fig : interpolated result of MAP_interpolator.m
        InterpolatedPowerMap(GP).fig : interpolated result of GP_interpolator.m
        InterpolatedPowerMap(GP).fig : interpolated result of LS_interpolator.m
        AggregatedMSE.fig : result of AggregatedMSETest.m
        RegionalMSE.fig : result of RegionalMSETest.m
        SNR.fig : result of SNRTest.m
