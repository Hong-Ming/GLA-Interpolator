function PlotSNRTest

% Parameter
%%%%%%%%%%%%%%%%%%%%%%%%% SimulationParameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seed_data = 3;                   % the random seed for generating receive signal power diagram
seed_sample = 1;                 % the random seed for generating sample points
num_of_measurement = 400;        % number of measurement
num_of_instance = 100;           % number of realizaion for GMRF graph learning
%%%%%%%%%%%%%%%%%%%%%%%%% GaussianProcessParameter %%%%%%%%%%%%%%%%%%%%%%%%
SNR = 14;                     % gaussian noise variance
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
%%%%%%%%%%%%%%%%%%%%%%%%% Other %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_of_instance_m = 1000;

% Define text for file name
SEED_SAMPLE = sprintf('seed=%d',seed_sample);
INSTANCE = sprintf('%d',num_of_instance_m);
PAR = sprintf('[dc,L0,eta,sigma_psi]=[%d,%d,%d,%d]',dc,Lo,eta,sigma_psi);
STATION = sprintf('[RNG,REZ,SRC]=[%d,%d,{%d,%d}]',Range,Resolution,Source(1),Source(2));
GRAPH = sprintf([LaplacianType Solver 'Degree=%d'],Degree);
MEASURE = sprintf('N=%d',num_of_measurement);

% Define file name and path
MAPFileName = ['MSE' INSTANCE '(MAP),' SEED_SAMPLE ',' PAR ',' STATION ',' GRAPH ',' MEASURE '.mat'];
MAPMISFileName = ['MSE' INSTANCE '(MAP_MIS),' SEED_SAMPLE ',' PAR ',' STATION ',' GRAPH ',' MEASURE '.mat'];
GPFileName = ['MSE' INSTANCE '(GP),' SEED_SAMPLE ',' PAR ',' STATION ',' GRAPH ',' MEASURE '.mat'];
GPMISFileName = ['MSE' INSTANCE '(GP_MIS),' SEED_SAMPLE ',' PAR ',' STATION ',' GRAPH ',' MEASURE '.mat'];
AnalyticalFileName = ['MSE' INSTANCE '(Analytical),' SEED_SAMPLE ',' PAR ',' STATION ',' GRAPH ',' MEASURE '.mat'];

MAPFilePath = fullfile('SNRTest/',MAPFileName);
MAPMISFilePath = fullfile('SNRTest/',MAPMISFileName);
GPFilePath = fullfile('SNRTest/',GPFileName);
GPMISFilePath = fullfile('SNRTest/',GPMISFileName);
AnalyticalFilePath = fullfile('SNRTest/',AnalyticalFileName);

LoadData = load(MAPFilePath);
MAP = mean(LoadData.MSE,2);
LoadData = load(MAPMISFilePath);
MAP_MIS = mean(LoadData.MSE,2);
LoadData = load(GPFilePath);
GP = mean(LoadData.MSE,2);
LoadData = load(GPMISFilePath);
GP_MIS = mean(LoadData.MSE,2);
LoadData = load(AnalyticalFilePath);
Analytical = mean(LoadData.MSE,2);
increment = LoadData.increment;

figure;
grid on
hold on
box on
plot(20:increment:5,GP,'-^','MarkerEdgeColor','red','MarkerSize',8,'LineWidth',2,'Color','r')
plot(20:increment:5,MAP,'-o','MarkerEdgeColor','blue','MarkerSize',8,'LineWidth',2,'Color','b')
plot(20:increment:5,MAP_MIS,'-s','MarkerEdgeColor','cyan','MarkerSize',8,'LineWidth',2,'Color','c')
plot(20:increment:5,GP_MIS,'-d','MarkerEdgeColor','magenta','MarkerSize',8,'LineWidth',2,'Color','m')
plot(20:increment:5,Analytical,'k','LineWidth',2)
% title('Aggregated MSE vs. \#Measure point','FontSize',18,'Interpreter','latex')
xlabel('SNR','FontSize',15)
ylabel('Aggregated MSE','FontSize',15)
legend('GP','GLA','GLA Mismatch','GP Mismatch','Analytical','FontSize',20)
ylim([5,20])
hold off


end