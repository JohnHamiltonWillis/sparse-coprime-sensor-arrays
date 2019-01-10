function [C, x, sz] = coarrayTotal(N1, M1, N2, M2, flag)%%%%%%%
%%%%Plots the coarray for a given pair of linear arrays--- Subarray 1 with
%%%% N1 sensors, M1 undersampling factor and Subarray 2 with N2 sensors,
%%%% M2 undersampling factor.  Assumes the first sensor is being shared
    close all;
    L1 = 0:M1:M1*(N1-1);
    L2 = 0:M2:M2*(N2-1);
    L = [L1 L2];
    L = unique(L);
    numSensors = length(L);
    D = [];
    for idx = 1:numSensors-1
        for jdx = idx+1:numSensors
            D = [D L(idx)-L(jdx)];
        end
    end
    D = [D -D];
    D = [zeros(1,numSensors) D];
    D = sort(D);
    u = unique(D);
    x = min(u):max(u);
    C = hist(D,x);
    if flag
        
        f = figure('WindowState','maximized');
        ax = axes('Parent', f, 'FontWeight', 'Bold', 'FontSize', 16,... 
        'Position',[0.267203513909224 0.11 0.496339677891654 0.815]);
        hold on;
        stem(ax, x,C,'LineWidth',2);
        xlabel('lags (k)','FontWeight','Bold','FontSize',16);
        ylabel('Number of occurrences','FontWeight','Bold','FontSize',16);
        xlim([min(u) max(u)]);
        grid on;
    end
    Cplus = C((length(C)-1)/2+1:end);
    sz = find(Cplus==0, 1, 'first') - 1;
%     xlim([-(sz-1) (sz-1)]);

end
