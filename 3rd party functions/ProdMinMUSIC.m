function ProdMinMUSIC(M,N,U1,U2,add1,add2,varargin) %%%%%%%%%%%%%%%%%

%%%The first two arguments are M and N. U1 and U2 are the undersampling
%%%factors. add1 and add2 are the numbers of sensors added to subarray1 and
%%%subarray2 respectively. The first argument in varargin is for the list
%%%of directions in terms of cos(theta). The second argument in varargin is
%%%the number of snapshots. The third argument in varargin is SNR in dB.
    m = length(varargin);
    OptArgs = {[0.1 -0.6 0.5 0.9 0.77 0.33 -0.29 0.6 -0.5 ] 5000 0};
    OptArgs(1:m) = varargin;
    [us, SampleSize, SNRdB] = OptArgs{:};
    lambda = 50;
    d = lambda/2;
    kx = 2*pi/lambda * us;
    numSources = length(us);%%%%Number of sources
    %%%calculate noise variance for signal power 1
    vars = ones(size(us));
    varn = 10^(-SNRdB/10);
    indexA = (0:(M+add1-1))*U1;
    indexB = (0:(N+add2-1))*U2;
    maxIndex = max(indexA(end),indexB(end));
    L = maxIndex+1;%%%It is the number of sensors a filled array with the same aperture would have
    %%%proper Gaussian noise samples
    z = length(unique([indexA indexB]));
    nc = sqrt(varn/2)*randn(L,SampleSize) + 1i*sqrt(varn/2)*randn(L,SampleSize);
    %%%steering vector for the filler array
    vc = zeros(L,length(us));
    
    s = zeros(length(us), SampleSize);
    xc = zeros(size(vc,1), SampleSize);%%%%This is the total data for the whole array
    for idx = 1:length(us)
        %%%Calculate the direction vector corresponding to each signal
        vc(:,idx) = exp(-1i*kx(idx)*(0:maxIndex)*d);
        %%%input signals: proper Gaussian
        s(idx, :) = (sqrt(vars(idx)/2)*randn(1,SampleSize) + 1i*sqrt(vars(idx)/2)*randn(1,SampleSize));
        xc = xc + vc(:,idx)*s(idx,:);
    end
    %%%Add noise to the signal
    xc = xc + nc;%%%Now xc will have a matrix of data, but it will be like the transpose of the real data obtained by the studends
    x1 = xc(indexA+1,:);%%%Extract the data for Subarray1
    x2 = xc(indexB+1,:);%%%Extract the data for Subarray2
    %%output of the array
    u = -1:0.001:1;
    yprod = zeros(size(u));
    ymin = zeros(size(u));
    for idx = 1:length(u)
        %%array weights
        w1 = (exp(1i*2*pi/lambda * u(idx) *indexA.'*d))/(M+add1);
        w2 = (exp(1i*2*pi/lambda * u(idx) *indexB.'*d))/(M+add2);
        tempa = w1'*x1;
        tempb = w2'*x2;
        ymin(idx) = sum(min(abs([tempa;tempb])))/SampleSize;
        yprod(idx) = sum(tempa.*conj(tempb))/SampleSize;
    end
    ymin = ymin/max(abs(ymin));
    yprod = yprod/max(abs(yprod));
    %%%Find the autocorrelation at various lags. The second and third
    %%%arguments in the function ifourierTrans are the lowest lag and the
    %%%highest lag. The lowest should always be zero. But the highest lag
    %%%could be varied and it would determine the autocorrelation matrix
    %%%dimension. The highest lag should not exceed L-1 though.
    rmin = ifourierTrans((ymin.').^2,0,L-1);
    rprod = ifourierTrans((yprod.').^2,0,L-1);
    
    Rmin = toeplitz(rmin.');    
    Rprod = toeplitz(rprod.');    
    %%%%Find Noise matrix for Rmin
    [eVecd, eVald] = eig(Rmin);
    eVald = diag(eVald);
    [~,sortindexd] = sort(eVald,'descend');
    eVecsortedd = eVecd(:,sortindexd);
    noiseBasisd = eVecsortedd(:,length(us)+1:end);
    RNoisemin = noiseBasisd*noiseBasisd';
    %%%%Find Noise matrix for Rprod
    [eVecd, eVald] = eig(Rprod);
    eVald = diag(eVald);
    [~,sortindexd] = sort(eVald,'descend');
    eVecsortedd = eVecd(:,sortindexd);
    noiseBasisd = eVecsortedd(:,length(us)+1:end);
    RNoiseprod = noiseBasisd*noiseBasisd';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    count = 1;
    Pmin = zeros(size(u));
    Pprod = zeros(size(u));
    for idx = u  
        v = exp(1i*pi * idx*(0:(L-1)).');
        Pmin(count) = 1/(v'*RNoisemin*v);
        Pprod(count) = 1/(v'*RNoiseprod*v);
        count = count + 1;
    end
    Pmin = Pmin/max(abs(Pmin));
    Pprod = Pprod/max(abs(Pprod));    
    close all;
    f = figure;
    set(gcf, 'Position', [1 41 1366 651]); % Maximize figure.     
    ax = axes('Parent', f, 'FontWeight', 'Bold', 'FontSize', 16,... 
    'Position',[0.267203513909224 0.11 0.496339677891654 0.815]);
    hold all;
    plot(ax, u, 20*log10(abs(Pmin)), 'LineWidth', 2, 'Color', 'Black');
    hold on;
    plot(ax, u, 20*log10(abs(Pprod)), '-','LineWidth',2, 'Color', [0 0.6 0]);
    grid on;
    hold on;
    xlabel('cos(\theta)=u', 'FontSize', 16, 'FontWeight', 'Bold');
    ylabel('Output, dB', 'FontSize', 16, 'FontWeight', 'Bold');
    title(['Number of Sensors used = ' num2str(z), ' and Number of Sources = ' num2str(numSources),], 'FontSize', 16, 'FontWeight', 'Bold');
    legendhandle = legend('Min-MUSIC','Prod-MUSIC');
    set(legendhandle,'AutoUpdate','off');
    xlim([-1 1]);
    lowerlimit = -25;
    ylim([lowerlimit 0]);
    %%%%%The following loop marks the actual source locations by creating a
    %%%%%line corresponding to each source location
    for idx = 1:length(us)
        hold on;
        plot([us(idx) us(idx)],[lowerlimit 0],'r:','LineWidth',2);
    end   
end
