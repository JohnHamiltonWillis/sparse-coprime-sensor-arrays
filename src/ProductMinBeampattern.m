function [Bmin, Bprod, u, B1, B2] = ProductMinBeampattern(M, N, U1,U2,add1,add2)
%%% M and N are the number of sensors in the basic subarrays
%%% add1 and add2 are the additional sensors in subarray1 and subarray2
%%% U1 and U2 are undersampling factors
    uDelta = 0.001;
    u = -1:uDelta:1;
    %%%%The first argument in BeampatternLinearArray is set to 0 so that
    %%%%that function doesn't plot graphs. The second parameter is the
    %%%%number of sensors in a subarray. The third parameter is the number
    %%%%of extra sensors in a subarray for extension. The fourth parameter
    %%%%is the intersensor spacing in terms of lambda. So, need to specify
    %%%%0.5 as well.
    B1 = BeampatternLinearArray(0, M+add1, 0, U1*0.5);
    B1 = B1/max(abs(B1));
    B2 = BeampatternLinearArray(0, N+add2, 0, U2*0.5);
    B2 = B2/max(abs(B2));

    
    Bmin = min(abs([B1;B2]));
    Bmin = 20*log10(abs(Bmin));

    Bprod = B1.*conj(B2);
    Bprod = 10*log10(abs(Bprod));

    
    f1 = figure;
    a1 = axes('Parent', f1, 'FontSize', 16, 'FontWeight', 'Bold', ...
        'Position',[0.267203513909224 0.11 0.496339677891654 0.815]);
    hold(a1, 'all');
    box(a1, 'on');
    grid(a1, 'on');
    plot(a1, u, 20*log10(abs(B1)), 'LineStyle', '-.', 'LineWidth', 4, 'Color', 'Blue');
    hold(a1, 'all');
    plot(a1, u, 20*log10(abs(B2)), 'LineStyle', '-.', 'LineWidth', 4, 'Color', 'Red');
    hold(a1, 'all');

    plot(a1, u, Bmin, 'LineWidth', 3, 'Color', 'Black');
    hold on;
    plot(a1, u, Bprod, '--','LineWidth', 3, 'Color', [0 0.6 0]);
    hold on;
    plot(a1, [-1 1],[-13 -13],'--','LineWidth',1,'Color','k');
    xlabel('u=cos(\theta)', 'FontSize', 16, 'FontWeight', 'Bold');
    ylabel('20log|B(u)|, dB', 'FontSize', 16, 'FontWeight', 'Bold');
    set(gcf, 'Position', [50 1 1317 689]); % Maximicoprimepairs.pdfze figure.         
    legend('Subarray 1','Subarray 2','Min', 'Product');
    xlim([-1 1]);

    ylim([-30 0]);
    
 end
