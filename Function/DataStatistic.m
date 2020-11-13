function [S, u_xs, var_xs] = DataStatistic(SimulationParameter,GaussianProcessParameter,StationPositionParameter)

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


% Define parameter library
SimulationParameter = {seed_data,seed_sample,num_of_measurement,num_of_instance};
GaussianProcessParameter = {SNR, dc, Lo, eta, sigma_psi};
StationPositionParameter = {Range, Resolution, Source(1), Source(2)};

% Other parameter
sigma_n = sqrt(sigma_psi^2/(10^(SNR/10)));
N = (Range/Resolution)+1;        % plonts per row or column
num_of_points = N^2;             % number of points

% Define grid and coordination
x_coor = vec(repmat(0:Resolution:Range,Range/Resolution+1,1));
y_coor = repmat((Range:-Resolution:0)',Range/Resolution+1,1);
x_s = [x_coor y_coor];
difference_xasis = repmat(x_s(:,1),1,num_of_points)-repmat(x_s(:,1)',num_of_points,1);
difference_yasis = repmat(x_s(:,2),1,num_of_points)-repmat(x_s(:,2)',num_of_points,1);

% Define text for file name
TITLE = sprintf('X%d',num_of_instance);
GAUSSIAN_PAR = sprintf('[SNR,dc,L0,eta,sigma_psi]=[%d,%d,%d,%d,%d]',SNR,dc,Lo,eta,sigma_psi);
SEED_SAMPLE = sprintf('seed=%d',seed_sample);
SAMPLE = sprintf('sample=%d',num_of_measurement);
STATION = sprintf('[RNG,REZ,SRC]=[%d,%d,{%d,%d}]',Range,Resolution,Source(1),Source(2));

% Define file name
DataStatisticFile = [TITLE ',' GAUSSIAN_PAR ',' STATION '.mat'];
SampleDataFile = [SEED_SAMPLE ',' SAMPLE ',' STATION '.mat'];

% Define file path
DataStatisticPath = fullfile('Data/GraphLearningData/',DataStatisticFile);
SampleDataPath = fullfile('Data/Sample/',SampleDataFile);

% load sample
LoadData = load(SampleDataPath);
measure = LoadData.measure;

% Construct data statistic
fprintf('1. Construct Data Statistic... '); tic;
W = fprintf('(    /%4d)',num_of_instance);
X = zeros(num_of_points,num_of_instance);
if exist(DataStatisticPath,'file') == 2
    LoadData = load(DataStatisticPath);
    X = LoadData.X;
else
    for i = 1:num_of_instance
        fprintf(repmat('\b',1,W))
        fprintf('(%4d/%4d)',i,num_of_instance)
        seed_data = i;

        % Generate the true channel field (without noise)
        rng(seed_data)
        [GroundTrue,~,u] = MakeGroundTrue(GaussianProcessParameter,StationPositionParameter);

        % Generate observation channel field (adding gaussian noise with variance sigma_n)
        rng(seed_data)
        Observation = MakeObservation(GroundTrue,sigma_n);
        % store the measure points data to data statistic matrix
        xs = Observation(:);
        X(:,i) = xs;
    end
    save(DataStatisticPath,'X')
end
fprintf(repmat('\b',1,W))
ElapsedTime = toc;
fprintf('Done (Elapsed Time : %gs)\n',ElapsedTime);

X = X(measure,:);
u_xs = mean(X,2);
var_xs = mean(var(X,0,2));
X = X - repmat(u_xs,1,num_of_instance);
S = (1/num_of_instance)*(X*X');

end