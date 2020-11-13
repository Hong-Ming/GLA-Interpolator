function SNRTest
clc
% Parameter
%%%%%%%%%%%%%%%%%%%%%%%%% SimulationParameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seed_data = 3;                   % the random seed for generating receive signal power diagram
seed_sample = 1;                 % the random seed for generating sample points
num_of_measurement = 400;        % number of measure points
num_of_instance = 100;          % number of realizaion for GMRF graph learning
%%%%%%%%%%%%%%%%%%%%%%%%% GaussianProcessParameter %%%%%%%%%%%%%%%%%%%%%%%%
dc = 70;                         % correlation distance
Lo = 10;                         % path-loss constant
eta = 3;                         % path-loss exponent
sigma_psi = 9;                   % shadowing variance
%%%%%%%%%%%%%%%%%%%%%%%%% StationPositionParameter %%%%%%%%%%%%%%%%%%%%%%%%
Range = 200;                     % the range of the region
Resolution = 4;                  % the resolution of received signal diagram
Source = [100 100];                % the position of the base station
%%%%%%%%%%%%%%%%%%%%%%%%% GraphLearningParameter %%%%%%%%%%%%%%%%%%%%%%%%%%
LaplacianType = 'GGL';           % type of Laplacian matrix for GMRF (GGL, DDGL, CGL)
Solver = 'BCD';                  % type of solver
Degree = 2;                      % degree of polynimail for curve fitting 
%%%%%%%%%%%%%%%%%%%%%%%%% MachineLearningParameter %%%%%%%%%%%%%%%%%%%%%%%%
InitializeSeed = 1;              % random seed for initialization
lr_sigma_psi = 0.01;             % learning rate for sigma_psi
lr_dc = 0.001;                   % learning rate for dc
MaxIteration = 100;              % maximum iteration in inner loop
MaxCycle = 100;                  % maximum number of cycle
%%%%%%%%%%%%%%%%%%%%%%%%% Other %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Type = 'MAP_MIS';                     % type of model (MAP, GP or Analytical)
R = 1000;                        % number of monta carlo simulation
increment = -1;
SNR = 20:increment:5;
MSE = zeros(length(SNR),R);
Remake = false;

% Define parameter library
SimulationParameter = {seed_data,seed_sample,num_of_measurement,num_of_instance};
GaussianProcessParameter = {SNR(1), dc, Lo, eta, sigma_psi};
StationPositionParameter = {Range, Resolution, Source(1), Source(2)};
GraphLearningParameter = {LaplacianType,Solver,Degree};
MachineLearningParameter = {InitializeSeed, lr_sigma_psi, lr_dc, MaxIteration, MaxCycle};

% Define text for file name
SEED_SAMPLE = sprintf('seed=%d',seed_sample);
INSTANCE = sprintf('%d',R);
PAR = sprintf('[dc,L0,eta,sigma_psi]=[%d,%d,%d,%d]',dc,Lo,eta,sigma_psi);
STATION = sprintf('[RNG,REZ,SRC]=[%d,%d,{%d,%d}]',Range,Resolution,Source(1),Source(2));
GRAPH = sprintf([LaplacianType Solver 'Degree=%d'],Degree);
MEASURE = sprintf('N=%d',num_of_measurement);

% Define file name and path
SNRFileName = ['MSE' INSTANCE '(' Type '),' SEED_SAMPLE ',' PAR ',' STATION ',' GRAPH ',' MEASURE '.mat'];
SNRFilePath = fullfile('MonteCarlo/SNRTest/',SNRFileName);

% make sample
[measure,~,~,~] = MakeSample(SimulationParameter,StationPositionParameter);

% Resume
if exist(SNRFilePath,'file') == 2 && ~Remake
    LoadData = load(SNRFilePath);
    r = LoadData.i;
    r = r + 1;
    MSE_old = LoadData.MSE;
    MSE(1:r-1,:) = MSE_old(1:r-1,:);
    if r > length(SNR)
        fprintf('Simulation Finished!\n')
    else
        fprintf('Resume from SNR = %d\n',SNR(r))
    end
else
    r = 1;
end

% Start NMSE test
if strcmp(Type,'Analytical')
    fprintf('Construct analytical variance...\n'); tic;
    num_of_points = ((Range/Resolution)+1)^2;
    [~,C,~] = MakeGroundTrue(GaussianProcessParameter,StationPositionParameter);
    for i = r:length(SNR)
        fprintf('SNR = %d\n',SNR(i))
        sigma_n = sqrt(sigma_psi^2/(10^(SNR(i)/10)));
        Sigma_RX = zeros(num_of_points,1);
        K = C + sigma_n^2 * eye(size(C));
        Ks = K(measure,measure);
        for j = 1:num_of_points
            ks = C(measure,j);
            Sigma_RX(j) = sigma_psi^2 - ks'*(Ks\ks);
        end
        MSE(i,:) = (1/num_of_points)*sum(Sigma_RX);
        save(SNRFilePath,'MSE','increment','i')
    end
else
    for i = r:length(SNR)
        GaussianProcessParameter{1} = SNR(i);
        for j = 1:R
            fprintf('# SNR = %d\n',SNR(i))
            SimulationParameter{1} = j+1000;
            if strcmp(Type,'MAP') || strcmp(Type,'MAP_MIS')
                [MSE(i,j),~,~,~] = GLA_Interpolator(SimulationParameter,GaussianProcessParameter,...
                                            StationPositionParameter,GraphLearningParameter);
            elseif strcmp(Type,'GP') || strcmp(Type,'GP_MIS')
                [MSE(i,j),~,~,~] = GP_Interpolator(SimulationParameter,GaussianProcessParameter,...
                                            StationPositionParameter,MachineLearningParameter);
            end
            clc;
        end
        save(SNRFilePath,'MSE','increment','i')
    end
end






end