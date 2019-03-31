function Subarray = CoprimeArray(M,N,U1,U2)
% Function to generate vector representation of the coprime array given Me,
% Ne, and the undersamping factors. Also returns the subarrays and the
% available lags and number of each lag. 
%%
% M = 2;
% N = 3;
% max_sensor = 64;

 max_sensor = max((M-1)*U1, (N-1)*U2)+1;
 
 
% Generate 0 1 representation of subarray 1 and 2
Subarray1 = zeros(1,max_sensor);
Subarray1((0:U1:(M-1)*U1)+1) = 1;
% Subarray1 = [1,Subarray1];
Subarray2 = zeros(1,max_sensor);
Subarray2((0:U2:(N-1)*U2)+1) = 1;
% Subarray2 = [1,Subarray2];

Sensor_placement = max(vertcat(Subarray1, Subarray2)); % Combine two subarrays

% Save the generated data to the Subarray struct
Subarray.array = Sensor_placement;
Subarray.sub1 = Subarray1;
Subarray.sub2 = Subarray2;

% Find the available lags given our coprime array
coarray = conv(Sensor_placement,fliplr(Sensor_placement));
Subarray.coarray = coarray;

end