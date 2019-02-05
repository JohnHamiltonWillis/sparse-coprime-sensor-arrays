function [coprime_array] = CoprimeLayoutOptimizer(sensor_layout,sensors)
% Function to use the best possible coprime pairs to reduce PSL height. 
% This function currently only generates coprime arrays starting at the
% first sensor in the array and only generates subarrays that have sensors
% available for full extension to the last sensor
cpairs = [];
coprime_spec = [];
for spacing = 1:63
    % Iterate through all the coprime pair spacings and generate pairs
    cpairs = GenerateCoprimePairs(2,64, spacing);
    for pair = 1:length(cpairs)
        for max_sensor = 2:64
            % Create the layout necessary given a comprime pair and number
            % of sensors
            subarray = CoprimeArray(cpairs{pair}(1),cpairs{pair}(2),max_sensor)
            % See if subarray for a coprime pair fits in available sensors
            if (sensor_layout & ~subarray)
                %pull PSL data
                temp = [2 3 13.25]; %temp = (pull data from table)
                % Store the highest PSL difference
                if temp(3) > coprime_spec(3)
                    coprime_spec = temp;
                end
            end
        end
    end
end
coprime_array = CoprimeArray(coprime_spec(1),coprime_spec(2),sensors);
end