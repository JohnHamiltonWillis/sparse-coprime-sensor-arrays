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