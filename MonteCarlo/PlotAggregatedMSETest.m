function PlotAggregatedMSETest

% Parameter
%%%%%%%%%%%%%%%%%%%%%%%%% SimulationParameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seed_sample = 1;                 % the random seed for generating sample points
num_of_instance = 100;           % number of realizaion for GMRF graph learning
%%%%%%%%%%%%%%%%%%%%%%%%% GaussianProcessParameter %%%%%%%%%%%%%%%%%%%%%%%%
SNR = 14;                        % gaussian noise variance
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
R = 1000;                        % number of monta carlo simulation

% Define text for file name
SEED_SAMPLE = sprintf('seed=%d',seed_sample);
INSTANCE = sprintf('%d',R);
PAR = sprintf('[SNR,dc,L0,eta,sigma_psi]=[%g,%d,%d,%d,%d]',SNR,dc,Lo,eta,sigma_psi);
STATION = sprintf('[RNG,REZ,SRC]=[%d,%d,{%d,%d}]',Range,Resolution,Source(1),Source(2));
GRAPH = sprintf([LaplacianType Solver 'Degree=%d'],Degree);

% Define file name and path
MAPFileName = ['MSE' INSTANCE '(MAP),' SEED_SAMPLE ',' PAR ',' STATION ',' GRAPH '.mat'];
MAPMISFileName = ['MSE' INSTANCE '(MAP_MIS),' SEED_SAMPLE ',' PAR ',' STATION ',' GRAPH '.mat'];
GPFileName = ['MSE' INSTANCE '(GP),' SEED_SAMPLE ',' PAR ',' STATION ',' GRAPH '.mat'];
GPMISFileName = ['MSE' INSTANCE '(GP_MIS),' SEED_SAMPLE ',' PAR ',' STATION ',' GRAPH '.mat'];
AnalyticalFileName = ['MSE' INSTANCE '(Analytical),' SEED_SAMPLE ',' PAR ',' STATION ',' GRAPH '.mat'];

MAPFilePath = fullfile('AggregatedMSETest/',MAPFileName);
MAPMISFilePath = fullfile('AggregatedMSETest/',MAPMISFileName);
GPFilePath = fullfile('AggregatedMSETest/',GPFileName);
GPMISFilePath = fullfile('AggregatedMSETest/',GPMISFileName);
AnalyticalFilePath = fullfile('AggregatedMSETest/',AnalyticalFileName);

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
plot(increment:increment:increment*length(GP),GP,'-^','MarkerEdgeColor','red','MarkerSize',8,'LineWidth',2,'Color','r')
plot(increment:increment:increment*length(MAP),MAP,'-o','MarkerEdgeColor','blue','MarkerSize',8,'LineWidth',2,'Color','b')
plot(increment:increment:increment*length(GP_MIS),GP_MIS,'-d','MarkerEdgeColor','magenta','MarkerSize',8,'LineWidth',2,'Color','m')
plot(increment:increment:increment*length(MAP_MIS),MAP_MIS,'-s','MarkerEdgeColor','cyan','MarkerSize',8,'LineWidth',2,'Color','c')
plot(increment:increment:increment*length(Analytical),Analytical,'k','LineWidth',2)
% title('Aggregated MSE vs. \#Measure point','FontSize',18,'Interpreter','latex')
xlabel('Number of measure points','FontSize',15)
ylabel('Aggregated MSE','FontSize',15)
legend('GP','GLA','GP Mismatch','GLA Mismatch','Analytical','FontSize',20)
xlim([100,600])
ylim([5,20])
hold off


end