close all;

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
periods = 4;
Array = CoprimeArray(N,M,4);
M_max = Array.M_N_max(1);
N_max = Array.M_N_max(2);
Max = Array.max+1;

sparseindices1 = find(Array.sub1);
sparseindices2 = find(Array.sub2);
%%%%To create x = (exp(1j*2*pi*f*t-1j*pi*cosd(-34)*indices)+ ...
%%%%exp(1j*2*pi*(-f)*t-1j*pi*cosd(-34)*indices)), make a matrix with the t
%%%%values and indices values first
times = (0:deltat:1)';
locations = 0:(Max-1);
[indices,t] = meshgrid(locations,times);
x = exp(1j*(2*pi*f*t-pi*cosd(angle)*indices));
clear indices t times;
%%%%To the matrix x above, we need to add white noise
totalData = x + sqrt(varn/2)*randn(size(x)) + 1i*sqrt(varn/2)*randn(size(x));

% sparseindices1 = (0:M:(M_max-1)) + 1;%%%We add 1 at the end because MATLAB...
                               %%%index starts at 1, not 0
% sparseindices2 = (0:N:(N_max-1)) + 1;%%%We add 1 at the end because MATLAB...
                               %%%index starts at 1, not 0

data1 = zeros(size(totalData));
data2 = zeros(size(totalData));
data1(:,sparseindices1) = totalData(:,sparseindices1);
data2(:,sparseindices2) = totalData(:,sparseindices2);

%% Process Data
calc_min = 1; % Set to one to calculate minimum processing output
calc_prod = 1; % ""                     product ""

% Min Processing
if calc_min == 1
    %%%%Find the product
    min_F = fft(data1,N*10,2) .* conj(fft(data2,N*10,2));
    %%%%Find the average
    min_mF = mean(abs(min_F));
    %%%%Normalize, fftshift, flip left to right and convert to dB
    min_mF = 10*log10(fliplr(fftshift(min_mF)/max(abs(min_mF))));
end

% Prod Processing
if calc_prod == 1;
    %%%%Find the minimum of the absolute values
    prod_F = min(abs(fft(data1,N*10,2)) , abs(fft(data2,N*10,2)));
    %%%%Find the average
    prod_mF = mean(abs(prod_F));
    %%%%Normalize, fftshift, flip left to right and convert to dB
    prod_mF = 20*log10(fliplr(fftshift(prod_mF)/max(abs(prod_mF))));
end

%% Plotting

plot_min = 1; % plot minimum processing output
plot_prod = 1; % ""  product ""

if calc_min == 1 && plot_min == 1
    w_min = linspace(-1,1,length(min_mF));
    figure_min = figure('WindowState','maximized');
    axes1 = axes('Parent',figure_min,...
        'Position',[0.13 0.11 0.775 0.331]);
    hold(axes1,'on');
    plot(w_min,min_mF,'Parent',axes1,'LineWidth',2);
    ylabel('Power dB','FontWeight','bold');
    xlabel('cos(\theta)','FontWeight','bold');
    title('Spatial Spectral Estimation','FontWeight','bold');
    xlim(axes1,[-1 1]);
    ylim(axes1,[-20 0]);
    box(axes1,'on');
    grid(axes1,'on');
    set(axes1,'FontSize',16,'FontWeight','bold');

    %%%%%Now evaluate the temporal spectrum
    F = fft(totalData,size(totalData,1)*10);
    mF = mean(abs(F),2);
    mF = fftshift(mF)/max(abs(mF));
    w = linspace(-1,1,length(mF));

    %%%%%%Plot the temporal spectrum
    axes2 = axes('Parent',figure_min,...
        'Position',[0.13 0.61 0.775 0.331]);
    hold(axes2,'on');
    plot(w,20*log10(mF),'Parent',axes2,'LineWidth',2);
    ylabel('Power dB','FontWeight','bold');
    xlabel('\times 0.5f_s, frequency (Hz)','FontWeight','bold');
    title('Temporal Spectral Estimation','FontWeight','bold');
    xlim(axes2,[-1 1]);
    ylim(axes2,[-40 0]);
    box(axes2,'on');
    grid(axes2,'on');
    % Set the remaining axes properties
    set(axes2,'FontSize',16,'FontWeight','bold');
end

if calc_prod == 1 && plot_prod == 1
    w_prod = linspace(-1,1,length(prod_mF));
    figure_prod = figure('WindowState','maximized');
    axes1 = axes('Parent',figure_prod,...
        'Position',[0.13 0.11 0.775 0.331]);
    hold(axes1,'on');
    plot(w_prod,prod_mF,'Parent',axes1,'LineWidth',2);
    ylabel('Power dB','FontWeight','bold');
    xlabel('cos(\theta)','FontWeight','bold');
    title('Spatial Spectral Estimation','FontWeight','bold');
    xlim(axes1,[-1 1]);
    ylim(axes1,[-20 0]);
    box(axes1,'on');
    grid(axes1,'on');
    set(axes1,'FontSize',16,'FontWeight','bold');

    %%%%%Now evaluate the temporal spectrum
    F = fft(totalData,size(totalData,1)*10);
    mF = mean(abs(F),2);
    mF = fftshift(mF)/max(abs(mF));
    w = linspace(-1,1,length(mF));

    %%%%%%Plot the temporal spectrum
    axes2 = axes('Parent',figure_prod,...
        'Position',[0.13 0.61 0.775 0.331]);
    hold(axes2,'on');
    plot(w,20*log10(mF),'Parent',axes2,'LineWidth',2);
    ylabel('Power dB','FontWeight','bold');
    xlabel('\times 0.5f_s, frequency (Hz)','FontWeight','bold');
    title('Temporal Spectral Estimation','FontWeight','bold');
    xlim(axes2,[-1 1]);
    ylim(axes2,[-40 0]);
    box(axes2,'on');
    grid(axes2,'on');
    % Set the remaining axes properties
    set(axes2,'FontSize',16,'FontWeight','bold');
end