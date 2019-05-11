function BeampatternsLinear
    M = 4; N = 5;e = 6;
    Me = e*M;
    Ne = e*N;
    ind1 = (0:(Me-1))*N;
    ind2 = (0:(Ne-1))*M;
    nind = unique([ind1, ind2]);
    u = -1:0.001:1;
    
        
    B1 = zeros(size(u));
    B2 = zeros(size(u));
    Bn = zeros(size(u));
    for idx = 1:length(u)
        v1 = exp(1i*pi*u(idx)*ind1')/Me;
        B1(idx) = sum(v1);
        v2 = exp(1i*pi*u(idx)*ind2')/Ne;
        B2(idx) = sum(v2);
        vn = exp(1i*pi*u(idx)*nind')/(length(nind));
        Bn(idx) = sum(vn);
    end
    close all;
        f = figure;
        set(gcf, 'Position', [50 1 1317 689]); % Maximize figure.     
        a = axes('Parent', f, 'FontWeight', 'Bold', 'FontSize', 16, ...
            'Position',[0.267203513909224 0.11 0.496339677891654 0.815]);
%         plot(a, u, 20*log10(abs(B1)), 'LineWidth', 3, 'Color', 'r');
%         hold on;
%         plot(a, u, 20*log10(abs(B2)), 'LineWidth', 3, 'Color', 'b');
%         hold on;
        plot(a, u, 10*log10(abs(B1.*conj(B2))), ':','LineWidth', 3, 'Color', 'k');
        hold on;
        plot(a, u, 20*log10(abs(Bn)),'-.', 'LineWidth', 3, 'Color', 'c');
        hold on;
        
        grid on;
        xlabel('u', 'FontSize', 16, 'FontWeight', 'Bold');
        ylabel('20log|B(u)|, dB', 'FontSize', 16, 'FontWeight', 'Bold');
         ylim([-30 0]);
        xlim([-1 1]);
end