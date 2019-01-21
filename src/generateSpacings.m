function Subarray = generateSpacings(M,N,periods)
% Function to generate vectors of the final sensor placements, as well as
% the placement for the subarrays. starts at sensor one since sensor zero
% is always there.
%%
% M = 3;
% N = 4;
% periods = 2;


M_N_max = [(periods*N-1)*M,(periods*M-1)*N]; % Calculate the max of each subarray
max_sensor = max(M_N_max); % Highest sensor number

% Generate 0 1 representation of subarray 1 and 2
Subarray1 = zeros(1,max_sensor);
Subarray1(M:M:M_N_max(1)) = 1;
Subarray1 = [1,Subarray1];
Subarray2 = zeros(1,max_sensor);
Subarray2(N:N:M_N_max(2)) = 1;
Subarray2 = [1,Subarray2];

Sensor_placement = max(vertcat(Subarray1, Subarray2)); % Combine two subarrays

% Save the generated data to the Subarray struct
Subarray.array = Sensor_placement;
Subarray.max = max_sensor;
Subarray.sub1 = Subarray1;
Subarray.sub2 = Subarray2;

% Find the available lags given our coprime array
lags = zeros(1,max_sensor+1); 
for i = 1:length(Sensor_placement) 
    if Sensor_placement(i) == 1 
        for k = i:length(Sensor_placement) 
            if Sensor_placement(k) == 1
                lags(k-i+1) = lags(k-i+1)+1; % add 1 to the number of lags at k-i lag
            end
        end
    end
end
lags = [0:max_sensor; lags];
Subarray.lags = lags;

end