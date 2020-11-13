function AggregatedMSETest
clc
% Parameter
%%%%%%%%%%%%%%%%%%%%%%%%% SimulationParameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seed_data = 3;                   % the random seed for generating receive signal power diagram
seed_sample = 1;                 % the random seed for generating sample points
num_of_instance = 100;           % number of realizaion for GMRF graph learning
%%%%%%%%%%%%%%%%%%%%%%%%% GaussianProcessParameter %%%%%%%%%%%%%%%%%%%%%%%%
SNR = 14;                        % signal to noise ratio
dc = 70;                         % correlation distance
Lo = 10;                         % path-loss constant
eta = 3;                         % path-loss exponent
sigma_psi = 9;                   % shadowing variance
%%%%%%%%%%%%%%%%%%%%%%%%% StationPositionParameter %%%%%%%%%%%%%%%%%%%%%%%%
Range = 200;                     % the range of the region
Resolution = 4;                  % the resolution of received signal diagram
Source = [100 100];              % the position of the base station
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
sigma_n = sqrt(sigma_psi^2/(10^(SNR/10)));
increment = 25;
measurement = 25:increment:600;
MSE = zeros(length(measurement),R);
Remake = false;


% Define parameter library
SimulationParameter = {seed_data,seed_sample,measurement(1),num_of_instance};
GaussianProcessParameter = {SNR, dc, Lo, eta, sigma_psi};
StationPositionParameter = {Range, Resolution, Source(1), Source(2)};
GraphLearningParameter = {LaplacianType,Solver,Degree};
MachineLearningParameter = {InitializeSeed, lr_sigma_psi, lr_dc, MaxIteration, MaxCycle};

% Define text for file name
SEED_SAMPLE = sprintf('seed=%d',seed_sample);
INSTANCE = sprintf('%d',R);
PAR = sprintf('[SNR,dc,L0,eta,sigma_psi]=[%g,%d,%d,%d,%d]',SNR,dc,Lo,eta,sigma_psi);
STATION = sprintf('[RNG,REZ,SRC]=[%d,%d,{%d,%d}]',Range,Resolution,Source(1),Source(2));
GRAPH = sprintf([LaplacianType Solver 'Degree=%d'],Degree);

% Define file name and path
AggregatedMSEFileName = ['MSE' INSTANCE '(' Type '),' SEED_SAMPLE ',' PAR ',' STATION ',' GRAPH '.mat'];
AggregatedMSEFilePath = fullfile('MonteCarlo/AggregatedMSETest/',AggregatedMSEFileName);

% Resume
if exist(AggregatedMSEFilePath,'file') == 2 && ~Remake
    LoadData = load(AggregatedMSEFilePath);
    r = LoadData.i;
    r = r + 1;
    MSE_old = LoadData.MSE;
    MSE(1:r-1,:) = MSE_old(1:r-1,:);
    if r > length(measurement)
        fprintf('Simulation Finished!\n')
    else
        fprintf('Resume from # Measurement = %d\n',measurement(r))
    end
else
    r = 1;
end

% Start NMSE test
if strcmp(Type,'Analytical')
    num_of_points = ((Range/Resolution)+1)^2;
    [~,C,~] = MakeGroundTrue(GaussianProcessParameter,StationPositionParameter);
    fprintf('Construct analytical variance...'); tic;
    for i = r:length(measurement)
        fprintf('# Measurement = %d\n',measurement(i))
        SimulationParameter{3} = measurement(i);
        [measure,~,~,~] = MakeSample(SimulationParameter,StationPositionParameter);
        Sigma_RX = zeros(num_of_points,1);
        K = C + sigma_n^2 * eye(size(C));
        Ks = K(measure,measure);
        
        for j = 1:num_of_points
            ks = C(measure,j);
            Sigma_RX(j) = sigma_psi^2 - ks'*(Ks\ks);
        end
        MSE(i,:) = (1/num_of_points)*sum(Sigma_RX);
        save(AggregatedMSEFilePath,'MSE','increment','i')
    end
    ElapsedTime = toc;
    fprintf(' Done (Elapsed Time : %gs)\n',ElapsedTime)
else
    for i = r:length(measurement)
        SimulationParameter{3} = measurement(i);
        for j = 1:R
            fprintf('# Measurement = %d\n',measurement(i))
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
        save(AggregatedMSEFilePath,'MSE','increment','i')
    end
end

end