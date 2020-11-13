function [GroundTrue,C,u] = MakeGroundTrue(GaussianProcessParameter,...
         StationPositionParameter)
     
% Parameter
dc = GaussianProcessParameter{2};
Lo = GaussianProcessParameter{3};
eta = GaussianProcessParameter{4};
sigma_psi = GaussianProcessParameter{5};
Range = StationPositionParameter{1};
Resolution = StationPositionParameter{2};
Source(1) = StationPositionParameter{3};
Source(2) = StationPositionParameter{4};
N = (Range/Resolution)+1;        % points per row or column
num_of_points = N^2;             % number of points
     
% Define grid and coordination
x_coor = repmat(0:Resolution:Range,Range/Resolution+1,1);
x_coor = x_coor(:);
y_coor = repmat((Range:-Resolution:0)',Range/Resolution+1,1);
grid_coor = [x_coor y_coor];
difference_xasis = repmat(grid_coor(:,1),1,num_of_points)-repmat(grid_coor(:,1)',num_of_points,1);
difference_yasis = repmat(grid_coor(:,2),1,num_of_points)-repmat(grid_coor(:,2)',num_of_points,1);
Distance = sqrt(difference_xasis.^2 + difference_yasis.^2);

% Make ground true
square = (grid_coor-repmat(Source,num_of_points,1)).*(grid_coor-repmat(Source,num_of_points,1));
u = Lo - 10*eta*log10(sqrt(square(:,1)+square(:,2)));
u(N*(Source(1)/Resolution)+N-(Source(2)/Resolution)) = 0;
C = sigma_psi^2 * exp(-Distance/dc);
p = mvnrnd(u,C);
GroundTrue = reshape(p,N,N);

end