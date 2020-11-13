addpath BCD_Algorithm
addpath Function


if exist('Data/GroundTrueData','dir') ~= 7
    mkdir Data/GroundTrueData
end
if exist('Data/Sample','dir') ~= 7
    mkdir Data/Sample
end
if exist('Data/GraphLearningData','dir') ~= 7
    mkdir Data/GraphLearningData
end
if exist('Data/GP','dir') ~= 7
    mkdir Data/GP
end


if exist('MonteCarlo/AggregatedMSETest','dir') ~= 7
    mkdir MonteCarlo/AggregatedMSETest
end
if exist('MonteCarlo/RegionalMSETest','dir') ~= 7
    mkdir MonteCarlo/RegionalMSETest
end
if exist('MonteCarlo/RegionalMSETest/Analytical','dir') ~= 7
    mkdir MonteCarlo/RegionalMSETest/Analytical
end
if exist('MonteCarlo/RegionalMSETest/Empirical','dir') ~= 7
    mkdir MonteCarlo/RegionalMSETest/Empirical
end
if exist('MonteCarlo/SNRTest','dir') ~= 7
    mkdir MonteCarlo/SNRTest
end