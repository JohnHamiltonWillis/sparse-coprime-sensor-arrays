function [Sensor_placement, Subarray1, Subarray2] = ...
    generateSpacings(M,N,periods)
% Function to generate vectors of the final sensor placements, as well as
% the placement for the subarrays. starts at sensor one since sensor zero
% is always there.
%%
% M = 3;
% N = 4;
% periods = 2;
M_N_max = [(periods*N-1)*M,(periods*M-1)*N];
max_sensor = max(M_N_max);
Subarray1 = zeros(1,max_sensor);
Subarray1(M:M:M_N_max(1)) = 1;
Subarray2 = zeros(1,max_sensor);
Subarray2(N:N:M_N_max(2)) = 1;
Sensor_placement = max(vertcat(Subarray1, Subarray2));
end