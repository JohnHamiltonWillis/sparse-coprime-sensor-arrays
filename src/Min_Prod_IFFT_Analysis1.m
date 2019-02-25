% close all;
% clear;
% clc;

% Generate Data
angle = -34;
f = 8923.8;%%%%Signal frequency in Hz
c = 343;%%%%%Signal speed in m/s
deltat = 1/44100;%%%%Temporal sampling interval in s
SNRdB = -10;%%%%SNR in dB
vars = 1;%%%%Signal variance
varn = vars*10^(-SNRdB/10);%%%%Noise variance

N = 4;
M = 5;
max_sensor = 80;
Array = CoprimeArray(N,M,max_sensor);

max_sensor = Array.max_sensor;

sparseindices1 = find(Array.sub1);
sparseindices2 = find(Array.sub2);
%%%%To create x = (exp(1j*2*pi*f*t-1j*pi*cosd(-34)*indices)+ ...
%%%%exp(1j*2*pi*(-f)*t-1j*pi*cosd(-34)*indices)), make a matrix with the t
%%%%values and indices values first
times = (0:deltat:1)';
locations = 0:(max_sensor-1);
[indices,t] = meshgrid(locations,times);
x = exp(1j*(2*pi*f*t-pi*cosd(angle)*indices));
clear indices t times;
%%%%To the matrix x above, we need to add white noise
totalData = x + sqrt(varn/2)*randn(size(x)) + 1i*sqrt(varn/2)*randn(size(x));

data1 = zeros(size(totalData));
data2 = zeros(size(totalData));
data1(:,sparseindices1) = totalData(:,sparseindices1);
data2(:,sparseindices2) = totalData(:,sparseindices2);

%% Process Data
calc_min = 1; % Set to one to calculate minimum processing output
calc_prod = 1; % ""                     product ""

% Min Processing
if calc_min
    %%%%Find the minimum of the absolute values
    min_F = min(abs(fft(data1,max_sensor*10,2)) , abs(fft(data2,max_sensor*10,2)));
    %%%%Find the average
    min_mF = mean(abs(min_F));
    %%%%Normalize, fftshift, flip left to right and convert to dB
    min_mF = 20*log10(fliplr(fftshift(min_mF)/max(abs(min_mF))));
end
 
% Prod Processing
if calc_prod
    %%%%Find the product
    prod_F = fft(data1,max_sensor*10,2) .* conj(fft(data2,max_sensor*10,2));
    %%%%Find the average
    prod_mF = mean(abs(prod_F));
    %%%%Normalize, fftshift, flip left to right and convert to dB
    prod_mF = 10*log10(fliplr(fftshift(prod_mF)/max(abs(prod_mF))));
end


%% Plotting

plot_min = 1; % plot minimum processing output
plot_prod = 1; % ""  product ""
overlay = 1; % if 1 overlay with preference (calc_min and calc_prod == 1)
if overlay && calc_min && calc_prod
    % Plot min
    w = linspace(-1,1,length(min_mF));
    figure_min = figure('WindowState','maximized');
    axes_min = axes('Parent',figure_min,...
        'Position',[0.13 0.11 0.775 0.331]);
    hold(axes_min,'on');
    plot(w,min_mF,'Parent',axes_min,'LineWidth',2,'Color','b');
    ylabel('Power dB','FontWeight','bold');
    xlabel('cos(\theta)','FontWeight','bold');
    title('Spatial Spectral Estimation - Overlay','FontWeight','bold');
    xlim(axes_min,[-1 1]);
    ylim(axes_min,[-20 0]);
    box(axes_min,'on');
    grid(axes_min,'on');
    set(axes_min,'FontSize',16,'FontWeight','bold');
    plot(w,prod_mF,'Parent', axes_min,'LineWidth',2,'Color','r');
    legend('Minimum','Product');
    
    %%%%%Now evaluate the temporal spectrum
    F = fft(totalData,size(totalData,1)*10);
    mF = mean(abs(F),2);
    mF = fftshift(mF)/max(abs(mF));
    w = linspace(-1,1,length(mF));

    %%%%%%Plot the temporal spectrum
    axes_2_min = axes('Parent',figure_min,...
        'Position',[0.13 0.61 0.775 0.331]);
    hold(axes_min,'on');
    plot(w,20*log10(mF),'Parent',axes_2_min,'LineWidth',2);
    ylabel('Power dB','FontWeight','bold');
    xlabel('\times 0.5f_s, frequency (Hz)','FontWeight','bold');
    title('Temporal Spectral Estimation','FontWeight','bold');
    xlim(axes_2_min,[-1 1]);
    ylim(axes_2_min,[-40 0]);
    box(axes_2_min,'on');
    grid(axes_2_min,'on');
    % Set the remaining axes properties
    set(axes_2_min,'FontSize',16,'FontWeight','bold');
    
    
end

% Plot min output
if ~overlay && calc_min && plot_min
    w = linspace(-1,1,length(min_mF));
    figure_min = figure('WindowState','maximized');
    axes_min = axes('Parent',figure_min,...
        'Position',[0.13 0.11 0.775 0.331]);
    hold(axes_min,'on');
    plot(w,min_mF,'Parent',axes_min,'LineWidth',2);
    ylabel('Power dB','FontWeight','bold');
    xlabel('cos(\theta)','FontWeight','bold');
    title('Spatial Spectral Estimation','FontWeight','bold');
    xlim(axes_min,[-1 1]);
    ylim(axes_min,[-20 0]);
    box(axes_min,'on');
    grid(axes_min,'on');
    set(axes_min,'FontSize',16,'FontWeight','bold');

    %%%%%Now evaluate the temporal spectrum
    F = fft(totalData,size(totalData,1)*10);
    mF = mean(abs(F),2);
    mF = fftshift(mF)/max(abs(mF));
    w = linspace(-1,1,length(mF));

    %%%%%%Plot the temporal spectrum
    axes_2_min = axes('Parent',figure_min,...
        'Position',[0.13 0.61 0.775 0.331]);
    hold(axes_min,'on');
    plot(w,20*log10(mF),'Parent',axes_2_min,'LineWidth',2);
    ylabel('Power dB','FontWeight','bold');
    xlabel('\times 0.5f_s, frequency (Hz)','FontWeight','bold');
    title('Temporal Spectral Estimation','FontWeight','bold');
    xlim(axes_2_min,[-1 1]);
    ylim(axes_2_min,[-40 0]);
    box(axes_2_min,'on');
    grid(axes_2_min,'on');
    % Set the remaining axes properties
    set(axes_2_min,'FontSize',16,'FontWeight','bold');
end

% Plot Prod output
if ~overlay && calc_prod && plot_prod
    w = linspace(-1,1,length(prod_mF));
    figure_prod = figure('WindowState','maximized');
    axes_prod = axes('Parent',figure_prod,...
        'Position',[0.13 0.11 0.775 0.331]);
    hold(axes_prod,'on');
    plot(w,prod_mF,'Parent',axes_prod,'LineWidth',2);
    ylabel('Power dB','FontWeight','bold');
    xlabel('cos(\theta)','FontWeight','bold');
    title('Spatial Spectral Estimation','FontWeight','bold');
    xlim(axes_prod,[-1 1]);
    ylim(axes_prod,[-20 0]);
    box(axes_prod,'on');
    grid(axes_prod,'on');
    set(axes_prod,'FontSize',16,'FontWeight','bold');

    %%%%%Now evaluate the temporal spectrum
    F = fft(totalData,size(totalData,1)*10);
    mF = mean(abs(F),2);
    mF = fftshift(mF)/max(abs(mF));
    w = linspace(-1,1,length(mF));

    %%%%%%Plot the temporal spectrum
    axes_2_prod = axes('Parent',figure_prod,...
        'Position',[0.13 0.61 0.775 0.331]);
    hold(axes_prod,'on');
    plot(w,20*log10(mF),'Parent',axes_2_prod,'LineWidth',2);
    ylabel('Power dB','FontWeight','bold');
    xlabel('\times 0.5f_s, frequency (Hz)','FontWeight','bold');
    title('Temporal Spectral Estimation','FontWeight','bold');
    xlim(axes_2_prod,[-1 1]);
    ylim(axes_2_prod,[-40 0]);
    box(axes_2_prod,'on');
    grid(axes_2_prod,'on');
    % Set the remaining axes properties
    set(axes_2_prod,'FontSize',16,'FontWeight','bold');
end

