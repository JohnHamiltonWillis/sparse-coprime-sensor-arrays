function [coprimes] = GenerateCoprimePairs(min,max,spacing)
% min is the lower bound in the range being searched. Max is the upper
% bound. Spacing sets the desired difference between the two coprimes. 

if spacing < 1
    spacing = 1;
end
if min <= spacing
    min = spacing+2;
    disp(['min too low, setting min to spacing +2 (' num2str(min) ')']);
end
if min < 3
    min = 3;
    disp('min too low, setting min to 3...');
end
if min - spacing == 1
    min = min + 1;
    disp(['min - spacing is equal to one, increasing min by one']);
end
if min > max
    error('Max must be greater than min');
end
%%
coprimes = cell(1, max);
for N = min:max
    M = N-spacing;
    factors_N = unique(factor(N));
        
    if ~ismember(M,factors_N)
        factors_M = unique(factor(M));
    else
        factors_M = M;
    end
    common_factors = intersect(factors_N, factors_M);
    if isempty(common_factors)
        coprimes{N} = [M,N];
    end

end
coprimes = coprimes(~cellfun('isempty',coprimes));
