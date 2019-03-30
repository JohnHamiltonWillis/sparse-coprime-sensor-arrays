%%%%This script finds the average of the MSEs over 1000 trials and
%%%%plots them
    snrRange = -15:15;
    product = zeros(size(snrRange));
    minimum = zeros(size(snrRange));
    direct = zeros(size(snrRange));
    full = zeros(size(snrRange));

    b = {'dataForSS50SNRp15Bias'};
    
    for count = 1:length(snrRange)
        load([b{1} num2str(snrRange(count)+15)]);
        product(count) = mean(p);
        minimum(count) = mean(m);
        direct(count) = mean(d);
        full(count) = mean(f);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    sdelta = 5;
    ssRange = 20:sdelta:200;
    product2 = zeros(size(ssRange));
    minimum2 = zeros(size(ssRange));
    direct2 = zeros(size(ssRange));
    full2 = zeros(size(ssRange));
    filename = {'dataForSNRp10SS'};
    ss = ssRange(1);
    for count = 1:length(ssRange)
        load ([filename{1} num2str(ss)]);
        ss = ss+sdelta;
        product2(count) = mean(p);
        minimum2(count) = mean(m);
        direct2(count) = mean(d);
        full2(count) = mean(f);
    end    
    
    
    close all;
    f1 = figure;
    a1 = axes('Parent', f1, 'FontSize', 16, 'FontWeight', 'Bold', ...
        'Position',[0.3375 0.587763289869609 0.366 0.337]);
    hold(a1, 'all');
    box(a1, 'on');
    grid(a1, 'on');
    plot(a1, snrRange, product, 'LineStyle', '-', 'LineWidth', 2, ...
        'Marker','square','Color', [0.6 0 0.6]);
    hold on;
    plot(a1, snrRange, minimum, 'LineStyle', '-.', 'LineWidth', 2, ...
        'Marker','diamond','Color', [0 0.6 0]);
    hold on;
    plot(a1, snrRange, direct, 'LineStyle', '--', 'LineWidth', 2, ...
        'Marker','square','Color', [0 0 0.6]);
    hold on;
    plot(a1, snrRange, full, 'LineStyle', ':', 'LineWidth', 2, ...
        'Marker','v','Color', [0 0 0]);
    hold on;

    xlabel('SNR, dB', 'FontSize', 16, 'FontWeight', 'Bold');
    ylabel('MSE', 'FontSize', 16, 'FontWeight', 'Bold');
    legend('Product','Minimum','Direct','Full');
    xlim([snrRange(1) snrRange(end)]);
    %set(a1,'yscale','log');

    a2 = axes('Parent', f1, 'FontSize', 16, 'FontWeight', 'Bold', ...
        'Position',[0.3375 0.140421263791374 0.366 0.337]);
    hold(a2, 'all');
    box(a2, 'on');
    grid(a2, 'on');
    plot(a2, ssRange, product2, 'LineStyle', '-', 'LineWidth', 2, ...
        'Marker','square','Color', [0.6 0 0.6]);
    hold on;
    plot(a2, ssRange, minimum2, 'LineStyle', '-.', 'LineWidth', 2, ...
        'Marker','diamond','Color', [0 0.6 0]);
    hold on;
    plot(a2, ssRange, direct2, 'LineStyle', '--', 'LineWidth', 2, ...
        'Marker','square','Color', [0 0 1]);
    hold on;
    plot(a2, ssRange, full2, 'LineStyle', ':', 'LineWidth', 2, ...
        'Marker','v','Color', [0 0 0]);
    hold on;
    
    xlabel('Snapshots', 'FontSize', 16, 'FontWeight', 'Bold');
    ylabel('MSE', 'FontSize', 16, 'FontWeight', 'Bold');
    set(gcf, 'Position', [50 1 1317 689]); 
    legend('Product','Minimum','Direct','Full');
    xlim([ssRange(1) ssRange(end)]);
    set(gcf,'WindowState','maximized');
    %set(a2,'yscale','log');
    
    
  