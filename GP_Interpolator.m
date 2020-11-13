function [MSE_Output,GroundTrue_Output,Observation_Output,Interpolation_Output]...
         = GP_Interpolator(SimulationParameter,GaussianProcessParameter,...
         StationPositionParameter,MachineLearningParameter)
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB version 
%     please use version R2019a or later
% Input arguments
%     1.SimulationParameter :{seed_data,seed_sample,num_of_measurement,num_of_instance}
%           1.seed_data : the random seed for generating receive signal power diagram
%           2.seed_sample : the random seed for generating sample points
%           3.num_of_measurement : number of measure points
%           4.num_of_instance : number of time instance for Gaussian process
%     2.GaussianProcessParameter : {SNR, dc, Lo, eta, sigma_psi}
%           1.SNR : signal to noise ratio
%           2.dc : correlation distance
%           3.Lo : path-loss constant
%           4.eta : path-loss exponent
%           5.sigma_psi : shadowing variance
%     3.StationPositionParameter : {Range, Resolution, Source_x, Source_y}
%           1.Range : the range of the region
%           2.Resolution : the resolution of the receive signal diagram
%           3.Source_x : the position of the base station (x coordinate)
%           4.Source_y : the position of the base station (y coordinate)
%     4.MachineLearningParameter : {InitializeSeed, lr_sigma_psi, lr_dc, MaxIteration, MaxCycle}
%           1.InitializeSeed : random seed for initialization
%           2.lr_sigma_psi : learning rate for sigma_psi
%           3.lr_dc : learning rate for dc
%           4.MaxIteration : maximum iteration in inner loop
%           5.MaxCycle : maximum number of cycle
% Output arguments
%     1.MSE_Output : aggregated mean square error
%     2.GroundTrue_Output : ground ture received signal power
%     3.Observation_Output : observed received signal power
%     4.Interpolation_Output : interpolated received signal power
% Options
%     1.Data Generation
%           1.RemakeData : always regenerate data and store it into database
%     2.Learning
%           1.RelearnData : always relearn data
%     3.Plot Figure
%           1.PlotMeasurePoints : plot sensor position on figures
%           2.PlotTruePowerMap : plot true power map
%           3.PlotObservedPowerMap : plot observed power map
%           4.PlotInterpolatedPowerMap : plot interpolated power map
% Usage
%     This is a polymorphic function, which works for any combination of
%     input and output. 
%     Example of usage : 
%           [Out1,Out2,Out3,Out4] = GP_Interpolator(In1,In2,In3,In4)
%           [Out1,~   ,Out3] = GP_Interpolator(In1,In2,[] ,In4)
%           [~   ,Out2] = GP_Interpolator(In1,[] ,In3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameter
%%%%%%%%%%%%%%%%%%%%%%%%% SimulationParameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seed_data = 1011;                % the random seed for generating receive signal power diagram
seed_sample = 1;                 % the random seed for generating sample points
num_of_measurement = 400;        % number of measure points
num_of_instance = 100;           % number of time instance for GMRF graph learning
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
%%%%%%%%%%%%%%%%%%%%%%%%% MachineLearningParameter %%%%%%%%%%%%%%%%%%%%%%%%
InitializeSeed = 1;              % random seed for initialization
lr_sigma_psi = 0.01;             % learning rate for sigma_psi
lr_dc = 0.001;                   % learning rate for dc
MaxIteration = 100;              % maximum iteration in inner loop
MaxCycle = 100;                  % maximum number of cycle

% Mode
%%%%%%%%%%%%%%%%%%%%%%%%% Data Generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RemakeData = false;
%%%%%%%%%%%%%%%%%%%%%%%%% Learning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RelearnData = false;
%%%%%%%%%%%%%%%%%%%%%%%%% Plot Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PlotMeasurePoints = false;
PlotTruePowerMap = true;
PlotObservedPowerMap = true;
PlotInterpolatedPowerMap = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Do Not Change Anything Below This Line %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Redefine parameter for different combination of input
if nargin >= 1
    if ~isempty(SimulationParameter)
        seed_data = SimulationParameter{1};
        seed_sample = SimulationParameter{2};
        num_of_measurement = SimulationParameter{3};
        num_of_instance = SimulationParameter{4};
    end
end
if nargin >=2
    if ~isempty(GaussianProcessParameter)
        SNR = GaussianProcessParameter{1};
        dc = GaussianProcessParameter{2};
        Lo = GaussianProcessParameter{3};
        eta = GaussianProcessParameter{4};
        sigma_psi = GaussianProcessParameter{5};
    end
end
if nargin >= 3
   if ~isempty(StationPositionParameter)
        Range = StationPositionParameter{1};
        Resolution = StationPositionParameter{2};
        Source(1) = StationPositionParameter{3};
        Source(2) = StationPositionParameter{4};
    end 
end
if nargin >= 4
    if ~isempty(MachineLearningParameter)
        InitializeSeed = MachineLearningParameter{1};
        lr_sigma_psi = MachineLearningParameter{2};
        lr_dc = MachineLearningParameter{3};
        MaxIteration = MachineLearningParameter{4};
        MaxCycle = MachineLearningParameter{5};
    end 
end

% Define parameter library
SimulationParameter = {seed_data,seed_sample,num_of_measurement,num_of_instance};
GaussianProcessParameter = {SNR, dc, Lo, eta, sigma_psi};
StationPositionParameter = {Range, Resolution, Source(1), Source(2)};

% Other parameter
sigma_n = sqrt(sigma_psi^2/(10^(SNR/10)));
N = (Range/Resolution)+1;        % points per row or column
num_of_points = N^2;             % number of points
BaseStation = N*(Source(1)/Resolution)+N-(Source(2)/Resolution);
BaseStation_xcoor = N-(Source(1)/Resolution);
BaseStation_ycoor = (Source(2)/Resolution)+1;
Distance_Noise = zeros(num_of_points,num_of_points);
rng(1)
Coor_Noise = 4*rand(num_of_points,2)-2;

% Robustness
if(num_of_measurement > num_of_points)
    error(['The number of sample : ' num2str(num_of_measurement) ...
    ' must be less than or equal to the number of points available for sampling : '...
    num2str(num_of_points) ])
end
if Source(1) > Range || Source(2) > Range
    error(['Source : (' num2str(Source(1)) ',' num2str(Source(2)) ...
    ') should be lie within ('  num2str(Range) ',' num2str(Range) ')'])
end
if mod(Range/Resolution,1) ~= 0
   error(['Range/Resolution : ' num2str(Range/Resolution) 'must be integer']) 
end

% Define text for file name
SEED_DATA = sprintf('seed=%d',seed_data);
SEED_SAMPLE = sprintf('seed=%d',seed_sample);
SAMPLE = sprintf('sample=%d',num_of_measurement);
GROUNDTRUE_PAR = sprintf('[dc,L0,eta,sigma_psi]=[%d,%d,%d,%d]',dc,Lo,eta,sigma_psi);
STATION = sprintf('[RNG,REZ,SRC]=[%d,%d,{%d,%d}]',Range,Resolution,Source(1),Source(2));

% Define file name
GroundTrueDataFile = [SEED_DATA ',' GROUNDTRUE_PAR ',' STATION '.mat'];
SampleDataFile = [SEED_SAMPLE ',' SAMPLE ',' STATION '.mat'];
GPDataFile = [SAMPLE '.mat'];

% Define file path
GroundTrueDataPath = fullfile('Data/GroundTrueData/',GroundTrueDataFile);
SampleDataPath = fullfile('Data/Sample/',SampleDataFile);
GPDataPath = fullfile('Data/GP/',GPDataFile);

% Define grid and coordination
x_coor = repmat(0:Resolution:Range,Range/Resolution+1,1);
x_coor = x_coor(:);
y_coor = repmat((Range:-Resolution:0)',Range/Resolution+1,1);
grid_coor = [x_coor y_coor] + Coor_Noise;
difference_xasis = repmat(grid_coor(:,1),1,num_of_points)-repmat(grid_coor(:,1)',num_of_points,1);
difference_yasis = repmat(grid_coor(:,2),1,num_of_points)-repmat(grid_coor(:,2)',num_of_points,1);
DistanceToMeasurement = sqrt(difference_xasis.^2 + difference_yasis.^2);
DistanceToStation = DistanceToMeasurement(:,BaseStation);

% Generate the true channel field (without noise)
rng(seed_data)
GroundTrue = zeros(N,N);
if exist(GroundTrueDataPath,'file') == 2 && ~RemakeData
    fprintf('Load Ground True Data : ~/');
    fprintf(GroundTrueDataPath);
    fprintf('\n');
    LoadData = load(GroundTrueDataPath);
    GroundTrue = LoadData.GroundTrue;
    C = LoadData.C;
    u = LoadData.u;
else
    fprintf('Make Ground True Data : ~/');
    fprintf(GroundTrueDataPath);
    fprintf('\n');
    [GroundTrue,C,u] = MakeGroundTrue(GaussianProcessParameter,StationPositionParameter);
    save(GroundTrueDataPath,'GroundTrue','C','u')
end

% Generate observation channel field (adding gaussian noise with variance sigma_n)
rng(seed_data)
Observation = zeros(N,N);
Observation = MakeObservation(GroundTrue,sigma_n);

% Ramdomly sample points
rng(seed_sample)
if exist(SampleDataPath,'file') == 2 && ~RemakeData
    fprintf('Load Sample Data : ~/');
    fprintf(SampleDataPath);
    fprintf('\n');
    LoadData = load(SampleDataPath);
    measure = LoadData.measure;
    measure_coor = LoadData.measure_coor;
else
    fprintf('Make Sample Data : ~/');
    fprintf(SampleDataPath);
    fprintf('\n');
    [measure,unknown,measure_coor,unknown_coor] = MakeSample(SimulationParameter,StationPositionParameter);
    save(SampleDataPath,'measure','unknown','measure_coor','unknown_coor')
end

[S, xs_bar, var_xs] = DataStatistic(SimulationParameter,GaussianProcessParameter,StationPositionParameter);
H = [ones(num_of_points,1), -10*log10(DistanceToStation)];
Hs = H(measure,:);


% Gradient Ascent
if exist(GPDataPath,'file') == 2 && ~RelearnData
    LoadData = load(GPDataPath);
    sigma_psi_hat = LoadData.sigma_psi_best;
    dc_hat = LoadData.dc_best;
    Lo_hat = LoadData.Lo_best;
    eta_hat = LoadData.eta_best;
    sigma_n_hat = sigma_n;
else
    rng(InitializeSeed)
    sigma_psi_square_old = var_xs-sigma_n^2;      % initial value for sigma_psi^2
    dc_old = rand;                                % initial value for dc
    Lo_old = rand;                                % initial value for Lo
    eta_old = rand;                               % initial value for eta
    theta_old = [Lo_old;eta_old]; 
    Threshold = 1;
    Cycle_Criteria = 1e-5;

    Cycle = 0;
    BestCycle = Threshold;
    sigma_psi_best = 0;
    dc_best = 0;
    Lo_best = 0;
    eta_best = 0;
    while Threshold > Cycle_Criteria && Cycle < MaxCycle
        Old = [sigma_psi_square_old,dc_old,theta_old(1),theta_old(2)];
        Cycle = Cycle + 1;
        fprintf('--------------------------------------------------------------\n')
        fprintf('Cycle = %3d\n',Cycle)
        
        % For sigma_psi
        Iteration = 0;
        WordCount = 0;
        fprintf('Iteration/sigma_psi = ');
        while Iteration < MaxIteration
            Ds = DistanceToMeasurement(measure,measure);
            Ks = sigma_psi_square_old*exp(-Ds/dc_old) + sigma_n^2*eye(num_of_measurement);

            dK = num_of_instance*(Ks\S/Ks - inv(Ks));
            dsigma_psi_square = exp(-Ds/dc_old);
            
            sigma_psi_square_new = sigma_psi_square_old + lr_sigma_psi*(1/2)*trace(dK*dsigma_psi_square);

            sigma_psi_square_old = sigma_psi_square_new;
            Iteration = Iteration + 1;
            fprintf(repmat('\b',1,WordCount))
            WordCount = fprintf('%3d/%4.2f',Iteration,sqrt(sigma_psi_square_new));
        end
        fprintf('\n')

        % For dc
        Iteration = 0;
        WordCount = 0;
        fprintf('Iteration/dc = ');
        while Iteration < MaxIteration
            Ds = DistanceToMeasurement(measure,measure);
            Ks = sigma_psi_square_old*exp(-Ds/dc_old) + sigma_n^2*eye(num_of_measurement);
            
            dK = num_of_instance*(Ks\S/Ks - inv(Ks));
            ddc = (sigma_psi_square_old/(dc_old^2))*Ds.*exp(-Ds/dc_old);

            dc_new = dc_old + lr_dc*(1/2)*trace(dK*ddc);
            dc_old = dc_new;
            Iteration = Iteration + 1;
            fprintf(repmat('\b',1,WordCount))
            WordCount = fprintf('%3d/%4.2f',Iteration,dc_new);
        end
        fprintf('\n')

        %For theta
        Ds = DistanceToMeasurement(measure,measure)/dc_old;
        Ks = sigma_psi_square_old*exp(-Ds) + sigma_n^2*eye(num_of_measurement);
        theta_new = (Hs'*(Ks\Hs))\Hs'*(Ks\xs_bar);
        theta_old = theta_new;
        fprintf('[Lo,eta] = [%4.2f,%4.2f]',theta_new(1),theta_new(2)')
        fprintf('\n')

        % Cycle Terminate Condition
        New = [sigma_psi_square_new,dc_new,theta_new(1),theta_new(2)];
        Threshold = (1/MaxIteration)*norm((New-Old),2);
        if Threshold < BestCycle
            BestCycle = Threshold;
            sigma_psi_best = sqrt(sigma_psi_square_new);
            dc_best = dc_new;
            Lo_best = theta_new(1);
            eta_best = theta_new(2);
            save(GPDataPath,'sigma_psi_best','dc_best','Lo_best','eta_best')
        end
        if BestCycle < 1e-2
           lr_sigma_psi = 0.001;
           lr_dc = 0.001;
        end
        if BestCycle < 1e-3
           lr_sigma_psi = 0.0005;
           lr_dc = 0.0005;
           MaxIteration = 50;
           MaxCycle = 200;
        end
        if BestCycle < 1e-4
           lr_sigma_psi = 0.0001;
           lr_dc = 0.0001;
           MaxIteration = 25;
           MaxCycle = 300;
        end
        fprintf('Error = %4.2e\n', Threshold)
    end
    sigma_psi_hat = sigma_psi_best;
    dc_hat = dc_best;
    Lo_hat = Lo_best;
    eta_hat = eta_best;
    sigma_n_hat = sigma_n;
end

C = sigma_psi_hat^2*exp(-DistanceToMeasurement/dc_hat);
Ks = C(measure,measure) + sigma_n_hat*eye(num_of_measurement);
u_hat = Lo_hat - 10*eta_hat*log10(DistanceToStation);
u_hat(BaseStation) = 0;
y = Observation(:) - u_hat;
ys = y(measure);

x_hat = zeros(num_of_points,1);
for i = 1:num_of_points
    k = C(measure,i);
    x_hat(i) = u_hat(i) + k'*(Ks\ys);
end
Interpolation = reshape(x_hat,N,N);
MSE = (1/num_of_points)*norm(Interpolation-GroundTrue,'fro')^2;


% Set the power at the station to zero
GroundTrue(BaseStation_xcoor,BaseStation_ycoor) = inf;
Observation(BaseStation_xcoor,BaseStation_ycoor) = inf;
Interpolation(BaseStation_xcoor,BaseStation_ycoor) = inf;

% Define output if the user specify output arguments
if nargout > 0
   MSE_Output = MSE; 
   GroundTrue_Output = GroundTrue;
   Observation_Output = Observation;
   Interpolation_Output = Interpolation;
end

% Set the power at the station to zero
GroundTrue(BaseStation_xcoor,BaseStation_ycoor) = 0;
Observation(BaseStation_xcoor,BaseStation_ycoor) = 0;
Interpolation(BaseStation_xcoor,BaseStation_ycoor) = 0;

GroundTrue(GroundTrue > 0) = 0;
Observation(Observation > 0) = 0;
Interpolation(Interpolation > 0) = 0;

% Set color bar limit
% cmax = max(max(max(Observation)),max(max(Interpolation)));
% cmin = min(min(min(Observation)),min(min(Interpolation)));
% cmax = 10*ceil(cmax/10);
% cmin = 10*floor(cmin/10);
cmax = 0;
cmin = -80;
increment = -8;
Level = cmax+increment:increment:cmin;

% Plot true receive signal power diagram
if PlotTruePowerMap
    figure;
    hold on
    contourf(0:Resolution:Range,Range:-Resolution:0,GroundTrue,Level)
    f1 = plot(Source(1),Source(2),'kx','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2.5);
    if PlotMeasurePoints
        for i = 1:num_of_measurement
           plot(measure_coor(i,1),measure_coor(i,2),'w+','MarkerSize',6,'LineWidth',1.5);
        end
        f2 = plot(NaN,NaN,'k+');
        legend([f1,f2],'Access Point','Measure Point');
    else
        legend(f1,'Measure Point');
    end
    xticks(0:50:200)
    yticks(0:50:200)
    xlim([0,Range])
    ylim([0,Range])
    xlabel('x (m)','FontSize',15)
    ylabel('y (m)','FontSize',15)
    title('True Power Map','FontSize',15)
    colorbar
    caxis([cmin cmax])
    hold off
end

% Plot observed receive sgnal power diagram
if PlotObservedPowerMap
    figure;
    hold on
    contourf(0:Resolution:Range,Range:-Resolution:0,Observation,Level)
    plot(Source(1),Source(2),'kx','MarkerFaceColor','k','MarkerSize',5,'LineWidth',2.5)
    f1 = plot(Source(1),Source(2),'kx','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2.5);
    for i = 1:num_of_measurement
       plot(measure_coor(i,1),measure_coor(i,2),'w+','MarkerSize',6,'LineWidth',1.5);
    end
    f2 = plot(NaN,NaN,'k+');
    legend([f1,f2],'Access Point','Measure Point');
    xlim([0,Range])
    ylim([0,Range])
    title('Observed Power Map','FontSize',15)
    xlabel({['$[d_c,L_0,\eta,\sigma_\psi]$ = ['...
    num2str(dc) ',' num2str(Lo) ',' num2str(eta) ',' num2str(sigma_psi) ...
    ']' '\ \ SNR = ' num2str(SNR) 'dB']},'FontSize',15,'Interpreter','latex')
    colorbar
    caxis([cmin cmax])
    hold off
end

% Plot interpolated receive sgnal power diagram
if PlotInterpolatedPowerMap
    figure;
    hold on
    contourf(0:Resolution:Range,Range:-Resolution:0,Interpolation,Level)
    plot(Source(1),Source(2),'kx','MarkerFaceColor','k','MarkerSize',6,'LineWidth',1.5)
    f1 = plot(Source(1),Source(2),'kx','MarkerFaceColor','k','MarkerSize',6,'LineWidth',2.5);
    if PlotMeasurePoints
        for i = 1:num_of_measurement
           plot(measure_coor(i,1),measure_coor(i,2),'w+','MarkerSize',6,'LineWidth',1.5);
        end
        f2 = plot(NaN,NaN,'k+');
        legend([f1,f2],'Access Point','Measure Point');
    else
        legend(f1,'Access Point');
    end
    xlim([0,Range])
    ylim([0,Range])
    xticks(0:50:200)
    yticks(0:50:200)
    xlabel('x (m)','FontSize',15)
    ylabel('y (m)','FontSize',15)
    title('Interpolated Power Map (GP)','FontSize',15)
    xlabel(['\#Measure Points = ' num2str(num_of_measurement) ...
            '\ \ MSE = ' num2str(MSE,'%7.4f')],'FontSize',15,'Interpreter','latex')
    colorbar
    caxis([cmin cmax])
    hold off
end

end