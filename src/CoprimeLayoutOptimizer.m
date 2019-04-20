function [coprime_array] = CoprimeLayoutOptimizer(sensor_layout)
%Function to use the best possible coprime pairs and number of periods 
%given a sparse array layout indicating which sensors have available data
%with a 1 or 0 for none. Current function is supported with data for up to
%64 sensors
if length(sensor_layout) >= 64
    error('sensor_layout too long. Use 64 or less sensors.');
end
cpairs = [];
for spacing = 1:63
    cpairs = GenerateCoprimePairs(2,64, spacing);
    for pair = 1:length(cpairs)
        CoprimeArray(cpairs{pair}(1),cpairs{pair}(2),64)
        %then do sensor_layout + CoprimeArray' RESUME PROGRESS HERE
    end
end