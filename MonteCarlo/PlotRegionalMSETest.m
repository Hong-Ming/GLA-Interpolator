function PlotRegionalMSETest

% Parameter
%%%%%%%%%%%%%%%%%%%%%%%%%SimulationParameter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seed_data = 3;                   % the random seed for generating receive signal power diagram
seed_sample = 1;                 % the random seed for generating sample points
num_of_measurement = 400;        % number of measure points
num_of_instance = 100;          % number of realizaion for GMRF graph learning
%%%%%%%%%%%%%%%%%%%%%%%%%GaussianProcessParameter%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR = 14;                        % sigmal to noise ratio
dc = 70;                         % correlation distance
Lo = 10;                         % path-loss constant
eta = 3;                         % path-loss exponent
sigma_psi = 9;                   % shadowing variance
%%%%%%%%%%%%%%%%%%%%%%%%%StationPositionParameter%%%%%%%%%%%%%%%%%%%%%%%%%%
Range = 200;                     % the range of the region
Resolution = 4;                  % the resolution of received signal diagram
Source = [100 100];                % the position of the base station
%%%%%%%%%%%%%%%%%%%%%%%%%Other%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = (Range/Resolution)+1;        % points per row or column
num_of_instance_f = 100;
Type = 'MAP';


% Define file name
SampleDataFile = sprintf('seed=%d,sample=%d,[RNG,REZ,SRC]=[%d,%d,{%d,%d}].mat'...
                         ,seed_sample,num_of_measurement,Range,Resolution,Source(1),Source(2));
VarianceFile = ...
sprintf(['Var%d(' Type '),R=%d[SNR,dc,L0,eta,sigma_psi]=[%g,%d,%d,%d,%d],[RNG,REZ,SRC]=[%d,%d,{%d,%d}].mat']...
        ,num_of_instance_f,num_of_measurement,SNR,dc,Lo,eta,sigma_psi,Range,Resolution,Source(1),Source(2));
    
% Define file path
AnalyticalMSEPath = fullfile('RegionalMSETest/Analytical',VarianceFile);
EmpiricalMSEPath = fullfile('RegionalMSETest/Empirical',VarianceFile);
% load data
LoadData = load(AnalyticalMSEPath);
AnalyticalMSE = LoadData.Var;
LoadData = load(EmpiricalMSEPath);
EmpiricalMSE = LoadData.Var;

% make sample
SimulationParameter = {seed_data,seed_sample,num_of_measurement,num_of_instance};
StationPositionParameter = {Range, Resolution, Source(1), Source(2)};
[~,~,measure_coor] = MakeSample(SimulationParameter,StationPositionParameter);

cmax = inf;
cmin = -inf;
% cmax = 26;
% cmin = 0;
Level = 0:3:24;
% Plot analytical MSE
figure;
hold on
contourf(0:Resolution:Range,Range:-Resolution:0,AnalyticalMSE,Level)
plot(Source(1),Source(2),'kx','MarkerFaceColor','k','MarkerSize',5,'LineWidth',2.5)
for i = 1:num_of_measurement
   plot(measure_coor(i,1),measure_coor(i,2),'W*') 
end
xlim([0,Range])
ylim([0,Range])
title('Analytical Regional MSE','FontSize',15)
xlabel({['\#Measure Points = ' num2str(num_of_measurement)] ...
    ['$[d_c,L_0,\eta,\sigma_\psi]$ = ['...
    num2str(dc) ',' num2str(Lo) ',' num2str(eta) ',' num2str(sigma_psi) ...
    ']' '  SNR = ' num2str(SNR) 'dB']},'FontSize',15,'Interpreter','latex')
colorbar
caxis([cmin cmax])
hold off

% set the MSE at the station to be the average of its adjacency points
EmpiricalMSE(N-(Source(1)/Resolution),(Source(2)/Resolution)+1)...
 = (1/4)*(EmpiricalMSE(N-(Source(1)/Resolution)+1,(Source(2)/Resolution)+1)...
        + EmpiricalMSE(N-(Source(1)/Resolution)-1,(Source(2)/Resolution)+1)...
        + EmpiricalMSE(N-(Source(1)/Resolution),(Source(2)/Resolution)+1+1)...
        + EmpiricalMSE(N-(Source(1)/Resolution),(Source(2)/Resolution)+1-1));

% Plot empirical MSE
figure;
hold on
contourf(0:Resolution:Range,Range:-Resolution:0,EmpiricalMSE,Level)
plot(Source(1),Source(2),'kx','MarkerFaceColor','k','MarkerSize',6,'LineWidth',1.5)
for i = 1:num_of_measurement
   plot(measure_coor(i,1),measure_coor(i,2),'W+','MarkerSize',6,'LineWidth',1.5) 
end
xticks(0:50:200)
yticks(0:50:200)
xlim([0,Range])
ylim([0,Range])
xlabel('x (m)','FontSize',15)
ylabel('y (m)','FontSize',15)
title(['Regional MSE ( ' Type ' )'],'FontSize',15)
xlabel({['\#Measure Points = ' num2str(num_of_measurement)...
    '  \#Instance = ' num2str(num_of_instance_f)],...
    ['$[d_c,L_0,\eta,\sigma_\psi]$ = ['...
    num2str(dc) ',' num2str(Lo) ',' num2str(eta) ',' num2str(sigma_psi) ...
    ']' '  SNR = ' num2str(SNR) 'dB']},'FontSize',15,'Interpreter','latex')
colorbar
caxis([cmin cmax])
hold off

end