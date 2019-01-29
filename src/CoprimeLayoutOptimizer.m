function [coprime_array] = CoprimeLayoutOptimizer(sensor_layout)
%Function to use the best possible coprime pairs and number of periods 
%given a sparse array layout indicating which sensors have available data
%with a 1 or 0 for none. Current function is supported with data for up to
%64 sensors
for spacing = 1:63
    cpairs(spacing) = GenerateCoprimePairs(0,64, spacing)
end