function [B, u] = BeampatternLinearArray(flag, N, varargin)
%%%If flag is 1, plot the graphs. 
%%%varargin: 1. shift = 0 for broadside steered, else enter u you want the
%%%array to steer to 2. intersensor spacing in terms of wavelength 3.
%%%weights: If you want to use matlab windows, use 'periodic' flag
%%%whenever available Example: w = hann(N, 'periodic').' and use w while
%%%calling the function 4. 5. xlimits
    m = length(varargin);
    OptArgs = {0 0.5  -1 1 1/N*ones(1,N)};
    OptArgs(1:m) = varargin;
    [shift, d, lu, uu, w] = OptArgs{:};
    n = (0:N-1)';
    u = lu:0.001:uu;
    B = zeros(size(u));
    for idx = 1:length(u)
        v = exp(1i*(n-(N-1)/2)*2*pi*(u(idx)-shift)*d);
%                 v = exp(1i*(n)*2*pi*(u(idx)-shift)*d);
        B(idx) = conj(w)*v;
    end
    % Normalize by the absolute value at u = 0
    B = B/sum(w);
    if flag
        f = figure;
        set(gcf, 'Position', [50 1 1317 689]); % Maximize figure.     
        a = axes('Parent', f, 'FontWeight', 'Bold', 'FontSize', 16, ...
            'Position',[0.267203513909224 0.11 0.496339677891654 0.815]);
        plot(a, u, (real(B)), 'LineWidth', 3, 'Color', 'Black');
        hold on;
        plot(a, u, (imag(B)), 'LineWidth', 3, 'Color', 'red');
%         title('Beam Pattern in U - Space', 'FontWeight', 'Bold', ...
%             'FontSize', 16);
        grid on;
        xlabel('u', 'FontSize', 16, 'FontWeight', 'Bold');
        ylabel('20log|B(u)|, dB', 'FontSize', 16, 'FontWeight', 'Bold');
%         ylim([-30 0]);
        xlim([-1 1]);
    end
end
