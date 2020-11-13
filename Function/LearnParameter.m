function [Sigma_s, x, u_xs] = LearnParameter(SimulationParameter,GaussianProcessParameter,...
         StationPositionParameter,GraphLearningParameter,Coor_Noise)

% s_quiet = cvx_quiet(true);

seed_data = SimulationParameter{1};
seed_sample = SimulationParameter{2};
num_of_measurement = SimulationParameter{3};
num_of_instance = SimulationParameter{4};

SNR = GaussianProcessParameter{1};
dc = GaussianProcessParameter{2};
Lo = GaussianProcessParameter{3};
eta = GaussianProcessParameter{4};
sigma_psi = GaussianProcessParameter{5};

Range = StationPositionParameter{1};
Resolution = StationPositionParameter{2};
Source(1) = StationPositionParameter{3};
Source(2) = StationPositionParameter{4};

LaplacianType = GraphLearningParameter{1};
Solver = GraphLearningParameter{2};
Degree = GraphLearningParameter{3};

Remake = true;

%%%%%%%%%%%%%%%% Do Not Change Anything Below This Line %%%%%%%%%%%%%%%%%%%

% Define parameter library
SimulationParameter = {seed_data,seed_sample,num_of_measurement,num_of_instance};
GaussianProcessParameter = {SNR, dc, Lo, eta, sigma_psi};
StationPositionParameter = {Range, Resolution, Source(1), Source(2)};
GraphLearningParameter = {LaplacianType,Solver,Degree};

% Other parameter
N = (Range/Resolution)+1;        % plonts per row or column
num_of_points = N^2;             % number of points

% Define text for file name
TITLE = sprintf(['Parameter%d(' LaplacianType ' ' Solver ')'],num_of_measurement);
GAUSSIAN_PAR = sprintf('[SNR,dc,L0,eta,sigma_psi]=[%d,%d,%d,%d,%d]',SNR,dc,Lo,eta,sigma_psi);
STATION = sprintf('[RNG,REZ,SRC]=[%d,%d,{%d,%d}]',Range,Resolution,Source(1),Source(2));
INSTANCE = sprintf('instance=%d',num_of_instance);
SEED_SAMPLE = sprintf('seed=%d',seed_sample);
SAMPLE = sprintf('sample=%d',num_of_measurement);
STATION = sprintf('[RNG,REZ,SRC]=[%d,%d,{%d,%d}]',Range,Resolution,Source(1),Source(2));


% Define file name
FitParameterFile = [TITLE ',' GAUSSIAN_PAR ',' STATION ',' INSTANCE '.mat'];
SampleDataFile = [SEED_SAMPLE ',' SAMPLE ',' STATION '.mat'];

% Define file path
FitParameterPath = fullfile('Data/GraphLearningData/',FitParameterFile);
SampleDataPath = fullfile('Data/Sample/',SampleDataFile);

LoadData = load(SampleDataPath);
measure = LoadData.measure;

% Define grid and coordination
x_coor = repmat(0:Resolution:Range,Range/Resolution+1,1);
x_coor = x_coor(:);
y_coor = repmat((Range:-Resolution:0)',Range/Resolution+1,1);
grid_coor = [x_coor y_coor]+Coor_Noise;
difference_xasis = repmat(grid_coor(:,1),1,num_of_points)-repmat(grid_coor(:,1)',num_of_points,1);
difference_yasis = repmat(grid_coor(:,2),1,num_of_points)-repmat(grid_coor(:,2)',num_of_points,1);
DistanceToMeasurement = sqrt(difference_xasis.^2 + difference_yasis.^2);

% GMRF graph learning
[~, Sigma_s, u_xs] = GMRF(SimulationParameter,GaussianProcessParameter,...
             StationPositionParameter,GraphLearningParameter);

% Start parameter learning
fprintf('3. Curving Fitting...'); tic;
D = DistanceToMeasurement;
Sigma_s_log = log(Sigma_s);
DSigma_s = -Sigma_s_log;
Dm = D(measure,measure);

if exist(FitParameterPath,'file') == 2 && ~Remake
    LoadData = load(FitParameterPath);
    x = LoadData.x;
else
    % Linear fit
    dm = zeros(length(Dm(:))-size(Dm,1),1);
    dSigma_s = zeros(length(DSigma_s(:))-size(DSigma_s,1),1);
    counter = 0;
    for i = 1:num_of_measurement
        for j = 1:num_of_measurement
            if i ~= j
               counter = counter + 1;
               dm(counter) = Dm(i,j);
               dSigma_s(counter) = DSigma_s(i,j);
            end
        end
    end
    B = ones(length(dm),Degree+1);
    for i = 1:Degree
        B(:,i+1) = dm.^i;
    end
    x = (B'*B)\B'*dSigma_s;
    save(FitParameterPath,'x')
end
ElapsedTime = toc;
fprintf(' Done (Elapsed Time : %gs)\n',ElapsedTime);

end