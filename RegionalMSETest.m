function RegionalMSETest

% Parameter
%%%%%%%%%%%%%%%%%%%%%%%%%SimulationParameter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seed_data = 3;                   % the random seed for generating receive signal power diagram
seed_sample = 1;                 % the random seed for generating sample points
num_of_measurement = 400;        % number of measure points
num_of_instance = 100;           % number of realizaion for GMRF graph learning
%%%%%%%%%%%%%%%%%%%%%%%%%GaussianProcessParameter%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR = 14;                        % signal to noise ratio
dc = 70;                         % correlation distance
Lo = 10;                         % path-loss constant
eta = 3;                         % path-loss exponent
sigma_psi = 9;                   % shadowing variance
%%%%%%%%%%%%%%%%%%%%%%%%%StationPositionParameter%%%%%%%%%%%%%%%%%%%%%%%%%%
Range = 200;                     % the range of the region
Resolution = 4;                  % the resolution of received signal diagram
Source = [100,100];              % the position of the base station
%%%%%%%%%%%%%%%%%%%%%%%%%Other%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = (Range/Resolution)+1;        % points per row or column
num_of_points = N^2;             % number of points
R = 1000;                        % number of monta carlo simulation

Type = 'MAP';                    % only for MAP
sigma_n = sqrt(sigma_psi^2/(10^(SNR/10)));

% Define parameter library
SimulationParameter = {seed_data,seed_sample,num_of_measurement,num_of_instance};
GaussianProcessParameter = {SNR, dc, Lo, eta, sigma_psi};
StationPositionParameter = {Range, Resolution, Source(1), Source(2)};

% Define file name
SampleDataFile = sprintf('seed=%d,sample=%d,[RNG,REZ,SRC]=[%d,%d,{%d,%d}].mat'...
                         ,seed_sample,num_of_measurement,Range,Resolution,Source(1),Source(2));
VarianceFile = ...
sprintf(['Var%d(' Type '),R=%d[SNR,dc,L0,eta,sigma_psi]=[%d,%d,%d,%d,%d],[RNG,REZ,SRC]=[%d,%d,{%d,%d}].mat']...
        ,R,num_of_measurement,SNR,dc,Lo,eta,sigma_psi,Range,Resolution,Source(1),Source(2));
    
% Define file path
SampleDataPath = fullfile('Data/Sample/',SampleDataFile);
AnalyticalVariancePath = fullfile('MonteCarlo/RegionalMSETest/Analytical',VarianceFile);
EmprircalMSEPath = fullfile('MonteCarlo/RegionalMSETest/Empirical',VarianceFile);

% Define grid and coordination
x_coor = repmat(0:Resolution:Range,Range/Resolution+1,1);
x_coor = x_coor(:);
y_coor = repmat((Range:-Resolution:0)',Range/Resolution+1,1);
grid_coor = [x_coor y_coor];
difference_xasis = repmat(grid_coor(:,1),1,num_of_points)-repmat(grid_coor(:,1)',num_of_points,1);
difference_yasis = repmat(grid_coor(:,2),1,num_of_points)-repmat(grid_coor(:,2)',num_of_points,1);
Distance = sqrt(difference_xasis.^2 + difference_yasis.^2);

% spatial autocovariance and mean
C = sigma_psi^2 * exp(-Distance/dc);
square = (grid_coor-repmat(Source,num_of_points,1)).*(grid_coor-repmat(Source,num_of_points,1));
u = Lo - 10*eta*log10(sqrt(square(:,1)+square(:,2)));
u(N*(Source(1)/Resolution)+N-(Source(2)/Resolution)) = 0;

% Ramdomly sample points
rng(seed_sample)
if exist(SampleDataPath,'file') == 2
    fprintf('Load Sample Data : ~/')
    fprintf(SampleDataPath)
    fprintf('\n')
    LoadData = load(SampleDataPath);
    measure = LoadData.measure;
else
    fprintf('Make Sample Data : ~/')
    fprintf(SampleDataPath)
    fprintf('\n')
    [measure,unknown,measure_coor,unknown_coor] = MakeSample(SimulationParameter,StationPositionParameter);
    save(SampleDataPath,'measure','unknown','measure_coor','unknown_coor')
end

% construct analytical variance
fprintf('Construct analytical variance...'); tic;
Var = zeros(N,N);
Sigma_RX = zeros(num_of_points,1);
K = zeros(num_of_measurement,num_of_measurement);
for i = 1:num_of_measurement
    for j = i:num_of_measurement
        K(i,j) = C(measure(i),measure(j));
        K(j,i) = C(measure(j),measure(i));
    end
end
K = K + sigma_n*eye(num_of_measurement);

for i = 1:num_of_points
    k = zeros(num_of_measurement,1);
    for j = 1:num_of_measurement
        k(j) = C(i,measure(j));
    end
    Sigma_RX(i) = sigma_psi^2 - k'*(K\k);
end
Var = reshape(Sigma_RX,N,N);
save(AnalyticalVariancePath,'Var')
ElapsedTime = toc;
fprintf(' Done (Elapsed Time : %gs)\n',ElapsedTime)

% Construct empirical MSE
fprintf('Construct empirical MSE...'); tic;
Var = zeros(N,N);
Mean = zeros(N,N);
X = zeros(N,N,R);
for i = 1:R
    SimulationParameter = {i,seed_sample,num_of_measurement,num_of_instance};
    [~,GroundTrue,~,Interpolation] = GLA_Interpolator(SimulationParameter,GaussianProcessParameter,...
         StationPositionParameter);
    Var = Var + (GroundTrue-Interpolation).^2;
    Mean = Mean + Interpolation;
    X(:,:,i) = Interpolation;
end
Var = Var / R;
save(EmprircalMSEPath,'Var')
ElapsedTime = toc;
fprintf(' Done (Elapsed Time : %gs)\n',ElapsedTime)

end