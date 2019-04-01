function [coprime_arrays] = CSAFinder(sensor_layout)
% sensor_layout is a vector of 1's and 0's where a 1 indicates that a 
% sensor is available at that position. array_length is the total number of
% sensor positions in the array. This function finds all the possible 
% coprime sensor array layouts that will fit in any given sparse linear 
% array.
coprime_arrays = [];
array_length = length(sensor_layout);
for spacing = 1:(array_length-1)
    % Iterate through all the coprime pair spacings and generate pairs
    cpairs = GenerateCoprimePairs(2,array_length, spacing);
    for pair = 1:length(cpairs)
        % Iterate through all the coprime pairs
        for max_sensor1 = 1:floor((array_length-1)/cpairs{pair}(1)+1)
            % Iterate through all the period extensions that fit in the
            % full sensor array for subarray1
            for max_sensor2 = 1:floor((array_length-1)/cpairs{pair}(2)+1)
                % Iterate through all the period extensions that fit in the
                % full sensor array for subarray2
                % Create the coprime layout necessary for a given coprime
                % pair and their subarrays extensions
                coprime_array = CoprimeArray(max_sensor1,max_sensor2,cpairs{pair}(1),cpairs{pair}(2));
                coprime_array = coprime_array.array;
                for shift = 1:(array_length-length(coprime_array))
                    % Shift coprime array from the beginning to the end of
                    % the full array
                    coprime_array = [0 coprime_array];
                    % Check if the proper sensors are available for coprime
                    % array
                    eligible = true;
                    for sensor = 1:length(coprime_array)             
                        if ((coprime_array(sensor)==1) && (sensor_layout(sensor)==0))
                            eligible = false;
                            break
                        end
                    end
                end
                % Store eligible coprime layout
                if eligible == true
                    coprime_arrays = [coprime_arrays; cpairs{pair}(1) cpairs{pair}(2) max_sensor1 max_sensor2];
                end
            end
        end
    end
end