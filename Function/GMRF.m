function [Theta, Sigma, u_xs] = GMRF(SimulationParameter,GaussianProcessParameter,...
                               StationPositionParameter,GraphLearningParameter)

%%%%%%%%%%%%%%%%%%%%%%%%% BCD parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = 0;
BCD_StoppingCriterion = 1e-9;
QP_StoppingCriterion = 1e-9;
MaxIteration = 300;
RegularizationType = 1;

%%%%%%%%%%%%%%%% Do Not Change Anything Below This Line %%%%%%%%%%%%%%%%%%%

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

% Define parameter library
SimulationParameter = {seed_data,seed_sample,num_of_measurement,num_of_instance};
GaussianProcessParameter = {SNR, dc, Lo, eta, sigma_psi};
StationPositionParameter = {Range, Resolution, Source(1), Source(2)};
GraphLearningParameter = {LaplacianType,Solver,Degree};

% Suppress output of CVX
if strcmp(Solver,'CVX')
    s_quiet = cvx_quiet(true);
end

% Define text for file name
TITLE = sprintf(['Graph%d(' LaplacianType ' ' Solver ')'],num_of_measurement);
GAUSSIAN_PAR = sprintf('[SNR,dc,L0,eta,sigma_psi]=[%d,%d,%d,%d,%d]',SNR,dc,Lo,eta,sigma_psi);
STATION = sprintf('[RNG,REZ,SRC]=[%d,%d,{%d,%d}]',Range,Resolution,Source(1),Source(2));
INSTANCE = sprintf('instance=%d',num_of_instance);

% Define file path
GraphFile = [TITLE ',' GAUSSIAN_PAR ',' STATION ',' INSTANCE '.mat'];
GraphFilePath = fullfile('Data/GraphLearningData/',GraphFile);

% Compute data statistic
[S, u_xs, ~] = DataStatistic(SimulationParameter,GaussianProcessParameter,StationPositionParameter);

% GMRF graph learning
fprintf('2. Solving Maximum Likelihood Estimation... '); tic;
if strcmp(Solver,'CVX')
    if exist(GraphFilePath,'file') == 2
        LoadData = load(GraphFilePath);
        Theta = LoadData.Theta;
        Sigma = LoadData.Sigma;
        status = LoadData.status;
    elseif strcmp(LaplacianType,'GGL')
        cvx_begin quiet
        variable Theta(num_of_measurement,num_of_measurement);
        minimize (trace(Theta*S) - log_det(Theta));
        subject to 
            Theta == semidefinite(num_of_measurement);
            for i = 1:num_of_measurement
                for j = 1:num_of_measurement
                    if i ~= j
                        Theta(i,j) <= 0;
                    end
                end
            end
        cvx_end
        Sigma = inv(Theta);
        status = cvx_status;
        save(GraphFilePath,'Theta','Sigma','status')
    elseif strcmp(LaplacianType,'DDGL')
        cvx_begin quiet
        variable Theta(num_of_measurement,num_of_measurement);
        minimize (trace(Theta*S) - log_det(Theta));
        subject to 
            Theta == semidefinite(num_of_measurement);
            sum(Theta,2) >= zeros(num_of_measurement,1);
            for i = 1:num_of_measurement
                for j = 1:num_of_measurement
                    if i ~= j
                        Theta(i,j) <= 0;
                    end
                end
            end
        cvx_end
        Sigma = inv(Theta);
        status = cvx_status;
        save(GraphFilePath,'Theta','Sigma','status')
    elseif strcmp(LaplacianType,'CGL')
        J = (1/num_of_measurement) * ones(num_of_measurement,num_of_measurement);
        cvx_begin quiet
        variable Theta(num_of_measurement,num_of_measurement);
        minimize (trace(Theta*(S+J)) - log_det(Theta+J));
        subject to 
            Theta == semidefinite(num_of_measurement);
            sum(Theta,2) >= -1e-7*ones(num_of_measurement,1);
            sum(Theta,2) <=  1e-7*ones(num_of_measurement,1);
            for i = 1:num_of_measurement
                for j = 1:num_of_measurement
                    if i ~= j
                        Theta(i,j) <= 0;
                    end
                end
            end
        cvx_end
        Sigma = inv(Theta);
        status = cvx_status;
        save(GraphFilePath,'Theta','Sigma','status')
    else
        error('Laplacian Type should be ''GGL'' or ''DDGL'' or ''CGL''')
    end
elseif strcmp(Solver,'BCD')
    Mask = ones(num_of_measurement,num_of_measurement);
    if exist(GraphFilePath,'file') == 2
        LoadData = load(GraphFilePath);
        Theta = LoadData.Theta;
        Sigma = LoadData.Sigma;
        status = LoadData.status;
    elseif strcmp(LaplacianType,'GGL')
        [Theta,Sigma,Info] = EstimateGGL(S,Mask,alpha,BCD_StoppingCriterion,QP_StoppingCriterion,...
                    MaxIteration,RegularizationType);
        if Info(end) <= BCD_StoppingCriterion
            status = 'Solved';
        else
            if Info(end) > BCD_StoppingCriterion * 1e3
                status = 'Fail';
            else
                status = 'Inaccurted Solved';
            end
        end    
        fprintf(repmat('\b',1,27))
        save(GraphFilePath,'Theta','Sigma','status')
    elseif strcmp(LaplacianType,'DDGL')
        [Theta,Sigma,Info] = EstimateDDGL(S,Mask,alpha,BCD_StoppingCriterion,QP_StoppingCriterion,...
                    MaxIteration,RegularizationType);
        if Info(end) <= BCD_StoppingCriterion
            status = 'Solved';
        else
            if Info(end) > BCD_StoppingCriterion * 1e3
                status = 'Fail';
            else
                status = 'Inaccurted Solved';
            end
        end    
        fprintf(repmat('\b',1,27))
        save(GraphFilePath,'Theta','Sigma','status')
    elseif strcmp(LaplacianType,'CGL')
        [Theta,Sigma,Info] = EstimateCGL(S,Mask,alpha,BCD_StoppingCriterion,QP_StoppingCriterion,...
                    MaxIteration,RegularizationType);
        if Info(end) <= BCD_StoppingCriterion
            status = 'Solved';
        else
            if Info(end) > BCD_StoppingCriterion * 1e3
                status = 'Fail';
            else
                status = 'Inaccurted Solved';
            end
        end    
        fprintf(repmat('\b',1,27))
        save(GraphFilePath,'Theta','Sigma','status')
    else
        error('Laplacian Type should be ''GGL'' or ''DDGL'' or ''CGL''')
    end
else
   error('Choose one solver') 
end
    
ElapsedTime = toc;
fprintf([status ' (Elapsed Time : %gs)\n'],ElapsedTime);

end