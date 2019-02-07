function Subarray = CoprimeArray(N,M,max_sensor)
% Function to generate vector representation of the coprime array given M,
% N, and the number of periods. Also returns the subarrays and the
% available lags and number of each lag. 
%%
% M = 2;
% N = 3;
% max_sensor = 64;



Subarray.max_sensor = max_sensor;

% Generate 0 1 representation of subarray 1 and 2
Subarray1 = zeros(1,max_sensor);
Subarray1((0:N:(max_sensor-1))+1) = 1;
% Subarray1 = [1,Subarray1];
Subarray2 = zeros(1,max_sensor);
Subarray2((0:M:(max_sensor-1))+1) = 1;
% Subarray2 = [1,Subarray2];


Sensor_placement = max(vertcat(Subarray1, Subarray2)); % Combine two subarrays

% Save the generated data to the Subarray struct
Subarray.array = Sensor_placement;
Subarray.sub1 = Subarray1;
Subarray.sub2 = Subarray2;

% Find the available lags given our coprime array
lags = zeros(1,max_sensor); 
for i = 1:length(Sensor_placement) 
    if Sensor_placement(i) == 1 
        for k = i:length(Sensor_placement) 
            if Sensor_placement(k) == 1
                lags(k-i+1) = lags(k-i+1)+1; % add 1 to the number of lags at k-i lag
            end
        end
    end
end
lags = [0:max_sensor-1; lags];
Subarray.lags = lags;

end