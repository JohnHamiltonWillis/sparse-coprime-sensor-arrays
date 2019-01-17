clear;close all;
%%%%%%% Some constant parameters
angle = -34;
N = 64;%%%%Number of sensors
f = 8923.8;%%%%Signal frequency in Hz
c = 343;%%%%%Signal speed in m/s
lambda = c/f;%%%%%Signal wavelength in m
deltat = 1/44100;%%%%Temporal sampling interval in s
SNRdB = -10;%%%%SNR in dB
vars = 1;%%%%Signal variance
varn = vars*10^(-SNRdB/10);%%%%Noise variance

%%%%To create x = (exp(1j*2*pi*f*t-1j*pi*cosd(-34)*indices)+ ...
%%%%exp(1j*2*pi*(-f)*t-1j*pi*cosd(-34)*indices)), make a matrix with the t
%%%%values and indices values first
times = (0:deltat:1)';
SampleSize = length(times);%%%%This is the number of samples per sensor
locations = 0:(N-1);
[indices,t] = meshgrid(locations,times);
x = cos(2*pi*f*t-pi*cosd(angle)*indices);
clear indices t times;
%%%%To the matrix x above, we need to add white noise
totalData = x + sqrt(varn/2)*randn(size(x)) + 1i*sqrt(varn/2)*randn(size(x));

%%%%%%%%%********************************************************************
%%%%%%%%%********************************************************************
%%%%%%%%%********************************************************************
%%%%%%%%%********************************************************************
%%%%%%%%%%%%%Uncomment the following line to use real data
% % load('/home/kayadhikari/Dropbox/research/AtTech/VDAMProject/Measurement_Data/4/2018-10-23_+23.mat');
% % totalData = totalData(1:40000,:);
%%%%%%%%%********************************************************************
%%%%%%%%%********************************************************************
%%%%%%%%%********************************************************************
%%%%%%%%%********************************************************************

%%%%The following four lines evaluate the spatial spectrum
F = fft(totalData,N*10,2);
mF = mean(abs(F));
mF = fliplr(fftshift(mF)/max(abs(mF)));
w = linspace(-1,1,length(mF));


%%%%%%%%%%%Plot the spatial spectrum
figure1 = figure('WindowState','maximized');
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.11 0.775 0.331]);
hold(axes1,'on');
plot(w,20*log10(mF),'Parent',axes1,'LineWidth',2);
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
axes2 = axes('Parent',figure1,...
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
