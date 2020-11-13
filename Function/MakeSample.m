function [measure,unknown,measure_coor,unknown_coor] = MakeSample(...
         SimulationParameter,StationPositionParameter)

% Parameter
seed_sample = SimulationParameter{2};
num_of_measurement = SimulationParameter{3};
Range = StationPositionParameter{1};
Resolution = StationPositionParameter{2};
Source(1) = StationPositionParameter{3};
Source(2) = StationPositionParameter{4};
N = (Range/Resolution)+1;
num_of_points = N^2;


rng(seed_sample)
% Make sample
measure = zeros(1,num_of_measurement);
unknown = zeros(1,num_of_points-num_of_measurement);
x_measure_coor = zeros(num_of_measurement,1);
y_measure_coor = zeros(num_of_measurement,1);
x_unknown_coor = zeros(num_of_points-num_of_measurement,1);
y_unknown_coor = zeros(num_of_points-num_of_measurement,1);

permuted_point = 1:num_of_points;
permuted_point(N*(Source(1)/Resolution)+N-(Source(2)/Resolution)) = num_of_points;
permuted_point(num_of_points) = N*(Source(1)/Resolution)+N-(Source(2)/Resolution);
if num_of_points ~= num_of_measurement
    for i = 1:num_of_measurement
    temp = permuted_point(i);
    index = randi([i num_of_points-1]);
    permuted_point(i) = permuted_point(index);
    permuted_point(index) = temp; 
    end
end
measure = permuted_point(1:num_of_measurement);
unknown = permuted_point(num_of_measurement+1:end);
measure = sort(measure);
unknown = sort(unknown);
for i = 1:num_of_measurement
x_measure_coor(i) = (floor((measure(i)-1) / N));
y_measure_coor(i) = N*(x_measure_coor(i)+1)-measure(i);
x_measure_coor(i) = x_measure_coor(i)*Resolution;
y_measure_coor(i) = y_measure_coor(i)*Resolution;
end
for i = 1:num_of_points-num_of_measurement
x_unknown_coor(i) = floor((unknown(i)-1) / N);
y_unknown_coor(i) = N*(x_unknown_coor(i)+1)-unknown(i);
x_unknown_coor(i) = x_unknown_coor(i)*Resolution;
y_unknown_coor(i) = y_unknown_coor(i)*Resolution;
end
measure_coor = [x_measure_coor, y_measure_coor];
unknown_coor = [x_unknown_coor, y_unknown_coor];

end