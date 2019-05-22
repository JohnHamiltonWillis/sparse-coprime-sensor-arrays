function [opt_coprime_array] = CoprimeLayoutOptimizer(sensor_layout)
% Function to use to determine the best possible coprime pairs to reduce
% PSL height. This function currently isn't hooked up to the data required 
% for comparisons. Once data is finished generating, this function will be
% updated.
opt_coprime_spec = [0 0 0];
for spacing = 1:63
    % Iterate through all the coprime pair spacings and generate pairs
    cpairs = GenerateCoprimePairs(2,64, spacing);
    for pair = 1:length(cpairs)
        % Iterate through all the coprime pairs
        for period1 = 1:floor(array_length/cpairs{pair}(1))
            % Iterate through all the period extensions that fit in the
            % full sensor array for subarray1
            max_sensor1 = period*cpairs{pair}(1);
            for period2 = 1:floor(array_length/cpairs{pair}(2))
                % Iterate through all the period extensions that fit in the
                % full sensor array for subarray2
                max_sensor2 = period*cpairs{pair}(2);
                % Create the coprime layout necessary for a given coprime
                % pair and their subarrays extensions
                coprime_array = CoprimeArray(max_sensor1,max_sensor2,cpairs{pair}(1),cpairs{pair}(2));
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
                % Store best coprime layout
                if eligible == true
                    % Pull PSL data
                    load(['C:\Users\Work\Downloads\2_100_' num2str(spacing) '.mat'], 'Z')
                    coprime_spec = Z(pair, period);
                    % Store the highest PSL difference
                    if coprime_spec(3) > opt_coprime_spec(3)
                        opt_coprime_spec = coprime_spec;
                    end
                end
            end
        end
    end
end
opt_coprime_array = CoprimeArray(opt_coprime_spec(1),opt_coprime_spec(2),64);
end

function [coprimes] = GenerateCoprimePairs(min,max,spacing)
% min is the lower bound in the range being searched. Max is the upper
% bound. Spacing sets the desired difference between the two coprimes. 

if spacing < 1
    spacing = 1;
end
if min < 2 % trivially don't want to test 1 since it's subarray is continuous
    min = 2;
    disp('Min too low, setting min to 2');
end
if min > max % Min must be greater than max
    error('Max must be greater than min');
end
%%
coprimes = cell(1, max); % Preallocate space 
for N = min+spacing:max % higher number n goes from min + spacing to max
    M = N-spacing; % lower number m
    factors_N = unique(factor(N)); % unique factors of N
        
    if ~ismember(M,factors_N) % if M is not a factor of N
        factors_M = unique(factor(M)); % find factors of M
    else
        factors_M = M; % else set M as the only factor of M
    end
    common_factors = intersect(factors_N, factors_M); % find intersection of factors of M and N
    if isempty(common_factors) % if there are no common factors
        coprimes{N} = [M,N]; % save pair as coprime
    end
end
coprimes = coprimes(~cellfun('isempty',coprimes)); % remove empty cell space
end

function Subarray = CoprimeArray(M,N,U1,U2)
% Function to generate vector representation of the coprime array given Me,
% Ne, and the undersamping factors. Also returns the subarrays and the
% available lags and number of each lag. 


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
Subarray.num_sensors = sum(Sensor_placement);

% Find the available lags given our coprime array
coarray = conv(Sensor_placement,fliplr(Sensor_placement));
Subarray.coarray = coarray;

end