function [opt_coprime_array] = CoprimeLayoutOptimizer(sensor_layout)
% Function to use the best possible coprime pairs to reduce PSL height. 
% This function currently only generates coprime arrays starting at the
% first sensor in the array and only generates subarrays that have sensors
% available for full extension to the last sensor
coprime_spec = [0 0 0];
for spacing = 1:63
    % Iterate through all the coprime pair spacings and generate pairs
    cpairs = GenerateCoprimePairs(2,64, spacing);
    for pair = 1:length(cpairs)
        for max_sensor = 3:64
            % Create the layout necessary given a comprime pair and number
            % of sensors
            coprime_array = CoprimeArray(cpairs{pair}(1),cpairs{pair}(2),max_sensor);
            coprime_array = coprime_array.array;
            % Shift the array across the full length of the array
            for shift = 0:(63-length(coprime_array))
                coprime_array = [0 coprime_array(1:(length(coprime_array)-1))];
            end
        end
    end
end
opt_coprime_array = CoprimeArray(coprime_spec(1),coprime_spec(2),64);