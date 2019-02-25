function [opt_coprime_array] = CoprimeLayoutOptimizer(sensor_layout)
% Function to use to determine the best possible coprime pairs to reduce
% PSL height. This function currently isn't hooked up to the data required 
% for comparisons. Once data is in a more readable format, this algorithm
% will be updated.
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
                % Check if the proper sensors are available for coprime
                % array
                eligible = true;
                for sensor = 1:length(coprime_array)             
                    if ((coprime_array(sensor)==1) && (sensor_layout(sensor)==0))
                        eligible = false;
                        break
                    end
                end
                % Store best coprime layout
                if eligible == true
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
end
opt_coprime_array = CoprimeArray(coprime_spec(1),coprime_spec(2),64);