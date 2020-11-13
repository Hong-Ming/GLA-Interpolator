# Graph Learning and Augmentation Based Interpolation of Signal Strength for Location-aware Communicationsr
#### Authors: [Hong-Ming, Chiu](https://github.com/Hong-Ming), [Carrson C. Fung](https://mcube.nctu.edu.tw/~cfung/), [Antonio Ortega](https://viterbi.usc.edu/directory/faculty/Ortega/Antonio)
#### [Link to Paper](https://www.eurasip.org/Proceedings/Eusipco/Eusipco2020/pdfs/0002150.pdf)
Code regarding GMRF Graph Learning referenced from: [https://github.com/STAC-USC/Graph_Learning](https://github.com/STAC-USC/Graph_Learning)
#### Cite
Please cite [our paper](https://www.eurasip.org/Proceedings/Eusipco/Eusipco2020/pdfs/0002150.pdf) if you use this code in your own work:

```
@inproceedings{chiu2020gla,
  title={Graph Learning and Augmentation Based Interpolation of Signal Strength for Location-aware Communicationsr},
  author={Hong-Ming Chiu and Carrson C. Fung and Antonio Ortega},
  booktitle={European Signal Processing Conference (EUSIPCO)},
  year={2020}
}
```

## Intorduction
This MATLAB program contains the code for the paper "Graph learning and augmentation based interpolation of signal strength for location-aware communications". The code for the interpolator using Gaussian Process [2] and Kriging [3] are also included to perform comparison and Monta Carlo simulation.

## Contents
1. In `/`:
    - **GLA_interpolator.m**: main function for Graph Learning and Augmentation based interpolator [1]
    - GP_interpolator.m : main function for Gaussian Process based interpolator [2]
    - LS_interpolator.m : main function for Kriging based interpolator [3]
    - AggregatedMSETest.m : Monta Carlo simulation for aggregated MSE
    - RegionalMSETest.m : Monta Carlo simulation for Regonal MSE
    - SNRTest.m : Monta Carlo simulation for SNR v.s.MSE
    - startup.m : add path and make directory, executed automatically at startup
    - cleanup.m : clean all *.mat file in Data/ directory
2. In `Function/`: 
    - MakeGroundTrue.m : generate ground true power map [4]
    - MakeObservation.m : generate observed power map
    - MakeSample.m : generate sample points
    - DataStatistic.m : generate data statistic
    - GMRF.m : GMRF graph learning [5]
    - LearnParameter.m : learn paramters for graph augmentation
3. In `BCD_Algorithm/` [5] : 
    - EstimateGGL.m : GGL graph learning
    - EstimateDDGL.m : DDGL graph learning
    - EstimateCGL.m : CGL graph learning
    - nonnegative_qp_solver.m : function for solving inner QP subproblem
    - update_sherman_morrison_diag.m : funciton for updating diagonal elements in DDGL graph learning
4. In `MonteCarlo/` : 
    - AggregatedMSETest/ : store outputs of AggregatedMSETest.m
    - RegionalMSETest/ :store outputs of RegionalMSETest.m
    - SNRTest/ :store outputs of SNRTest.m
    - PlotAggregatedMSETest.m : plot results of AggregatedMSETest.m
    - PlotRegionalMSETest.m : plot results of RegionalMSETest.m
    - PlotSNRTest.m : plot results of SNRTest.m
5. In `Data/` :
    - GP/ : store results of GP_interpolator.m
    - GraphLearningData/ : store results of DataStatistic.m, GMRF.m and LearnParameter.m
    - GroundTrueData/ : store results of MakeGroundTrue.m
    - Sample/ : store results of MakeSample.m
6. In `Figure/` :
    - TruePowerMap.fig : true power map generated from MakeGroundTrue.m
    - ObservedPowerMap.fig : noisy observation of power map plotted in TruePowerMap.fig
    - InterpolatedPowerMap(MAP).fig : interpolated result of MAP_interpolator.m
    - InterpolatedPowerMap(GP).fig : interpolated result of GP_interpolator.m
    - InterpolatedPowerMap(GP).fig : interpolated result of LS_interpolator.m
    - AggregatedMSE.fig : result of AggregatedMSETest.m
    - RegionalMSE.fig : result of RegionalMSETest.m
    - SNR.fig : result of SNRTest.m
## SETUP/INSTALLATION
1. Require CVX [6] installed
2. Search paths and required directories will be create automatically as
   startup.m will automatically executed at startup.
3. If your startup directory is not GLA/ or there are errors occurs due to 
   'Function is not found on MATLAB searth path' or 'No such file or directory',
   please execute startup.m manually.
        
## RERERENCE
[1] H.M. Chiu, C.C. Fung, and A. Ortega, "Graph learning and augmentation based interpolation of signal strength for location-aware communications," 2020.
    
[2] Ferris Brian, Hähnel Dirk and Fox Dieter, "Gaussian Processes for Signal Strength-Based Location Estimation," 2006.
    
[3] A. Serrano, B. Girault, and A. Ortega, "Geostatistical Data Interpolation using Graph Signal Spectral Prior," 2019.
    
[4] R. Di Taranto et al., "Location-aware communicationsfor 5g networks," IEEE Signal Processing Magazine, vol. 31, no. 6, pp. 102?112, 2014.
    
[5] H.E. Egilmez, E. Pavez, and A. Ortega, "Graph learn-ing from data under laplacian and structural constraints," IEEE Journal of Selected Topics in Signal processing, vol. 11, no. 6, pp. 825?841, 2017.
    
[6] Michael Grant and Stephen Boyd. CVX: Matlab software for disciplined convex programming, version 2.0 beta. http://cvxr.com/cvx, September 2013.
    
## AUTOR/CONTACT INFO
Name  : Hong-Ming, Chiu

Email : hongmingchiu2017@gmail.com

