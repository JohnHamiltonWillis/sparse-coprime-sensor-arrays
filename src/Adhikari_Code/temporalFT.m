function temporalFT(varargin)
%%%%Uncomment and evaluate the statements at the end
    cla;
    m = length(varargin);
    OptArgs = {57};
    OptArgs(1:m) = varargin;   
    fs = 44100;
    [sensor] = OptArgs{:};
    load('/home/kayadhikari/Dropbox/research/AtTech/VDAMProject/Measurement_Data/3/2018-10-23_+23.mat');
    rows = 500;
    y = totalData(1:225000,sensor);
    clear totalData;
    x = reshape(y,rows,450);
    F = fft(x,rows*10);
    mF = mean(abs(F),2);
    clear F;
    mF = mF/max(abs(mF));
    w = linspace(-1,1,length(mF));
    plot(w*fs/2,20*log10(abs(fftshift(mF))),'LineWidth',2,'Color','k');
    ylim([-60 0]);
    xlim([-fs/2 fs/2]);
    set(gca,'FontWeight','Bold','FontSize',16);
    set(gcf,'WindowState','maximized');
    grid on;
    hold on;
    title(['Sensor Number: ',num2str(sensor)]);
    
% % close all;
% % for idx = 1:64
% % temporalFT(idx);pause;
% % end    