function directionEstimatesVersion5
        SNRdB = 10;
        SampleSize = 100;
        plotfigure = 1;
        M = 4; N = 6;
        U1 = 3; U2 = 2;
        us = cosd(randi(181,[1 2])-1);
        numSources = length(us);
        lambda = 50;    d = lambda/2;    kx = 2*pi/lambda * us;
        %%%calculate noise variance for signal power 1
        vars = ones(1,numSources);
        varn = vars(1)*10^(-SNRdB/10);
        s = zeros(numSources,SampleSize);
        %%input signals: proper GaussianK0478768904
        for idx = 1:numSources
            s(idx,:) = (sqrt(vars(idx)/2)*randn(1,SampleSize) + 1i*sqrt(vars(idx)/2)*randn(1,SampleSize));
        end
        ApertureEnd = max([(M-1)*U1 (N-1)*U2]);%%%%The array starts at 0 and ends at 63
        
        %%steering vector
        v = zeros(ApertureEnd+1,1);
        %%%%%The vector will have the data for all sensors
        x = zeros(ApertureEnd+1,SampleSize);
        for idx = 1:numSources
            v(:,idx) = exp(1i*kx(idx)*(0:ApertureEnd)*d).';
            x = x + v(:,idx)*s(idx,:);
        end
        %%%Add proper Gaussian noise samples
        x = x + sqrt(varn/2)*randn(ApertureEnd+1,SampleSize) + 1i*sqrt(varn/2)*randn(ApertureEnd+1,SampleSize);

        indexa = (0:U1:(M-1)*U1).';    indexb = (0:U2:(N-1)*U2).';
        indexunion = unique([indexa' indexb']);
        %%%xa will be subarray 1 data. It will have zero values where
        %%%sensors are not present.
        xa = zeros(max(indexa)+1,SampleSize);
        %%%xb will be subarray 2 data. It will have zero values where
        %%%sensors are not present
        xb = zeros(max(indexb)+1,SampleSize);%%%This will be subarray 2 data
        xa(indexa+1,:) = x(indexa+1,:);%%%xa takes the data from x, but only where Subarray 1 has sensors
        xb(indexb+1,:) = x(indexb+1,:);%%%xb takes the data from x, but only where Subarray 2 has sensors
        xtotal(indexunion+1,:) = x(indexunion+1,:);%%%%this is the union of Subarray 1 and Subarray 2. This data will be
                                                   %%%%used in direct MUSIC
                                                   

        sensorindicator(indexunion+1) = 1;%%%%This will make a vector called sensorindicator that will have
                                          %%%%ones where Subarray1 or Subarray 2 have sensors and
                                          %%%%zeros where both Subarray1 and Subarray 2 don't have sensors
        coarray = conv(sensorindicator,fliplr(sensorindicator));%%%%coarray tells you how many different 
                                                                %%%%sensor  pairs there are for each lag                                                         %%%%
        lags = -max(indexunion):max(indexunion);
        %%%%We can remove the negative half of the coarray and lags because
        %%%%the information they have is redundant
        temp = (length(coarray)-1)/2;%%%%This is the number of negative elements that can be removed
        coarray(1:temp) = [];%%%%Removes the negative half of coarray
        lags(1:temp) = [];%%%%Removes the negative half of lags
        %%%%The coarray could have holes. Keep the longest hole-free and
        %%%%the associated lags.
        firstzeroindex = find(coarray==0,1);%%%finds the index of the first zero elementconv(a,b)
        %%%%if the firstzeroindex is not empty, remove the elements from
        %%%%coarray and lag starting at the first zero element
        if ~isempty(firstzeroindex)
            coarray(firstzeroindex:end) = [];
            lags(firstzeroindex:end) = [];
        end
        r = zeros(length(coarray),SampleSize);%%%%covariance estimates

        for kdx = 1:SampleSize
            dataset = xtotal(:,kdx);
            %%%%The convolution operation can actually be used to find
            %%%%autocorrelation as shown below, for each set of samples
            tempR = conv(dataset.',fliplr(conj(dataset.')));
            tempR(1:temp) = [];
            if ~isempty(firstzeroindex)
                tempR(firstzeroindex:end) = [];
            end
            r(:,kdx) = (tempR.')./(coarray.');
        end
        rnew = zeros(11,SampleSize);        
        for kdx = 1:SampleSize
            dataset = xtotal(:,kdx);
            rnew(1,kdx) = xtotal(1,kdx)*conj(xtotal(1,kdx));
            rnew(2,kdx) = xtotal(3,kdx)*conj(xtotal(4,kdx));
            rnew(3,kdx) = xtotal(3,kdx)*conj(xtotal(5,kdx));
            rnew(4,kdx) = xtotal(1,kdx)*conj(xtotal(4,kdx));
            rnew(5,kdx) = xtotal(1,kdx)*conj(xtotal(5,kdx));
            rnew(6,kdx) = xtotal(5,kdx)*conj(xtotal(10,kdx));
            rnew(7,kdx) = xtotal(1,kdx)*conj(xtotal(7,kdx));
            rnew(8,kdx) = xtotal(3,kdx)*conj(xtotal(10,kdx));
            rnew(9,kdx) = xtotal(1,kdx)*conj(xtotal(9,kdx));
            rnew(10,kdx) = xtotal(1,kdx)*conj(xtotal(10,kdx));
            rnew(11,kdx) = xtotal(1,kdx)*conj(xtotal(11,kdx));
            tempR = conv(dataset.',fliplr(conj(dataset.')));
            tempR(1:temp) = [];
            if ~isempty(firstzeroindex)
                tempR(firstzeroindex:end) = [];
            end
            r(:,kdx) = (tempR.')./(coarray.');
        end
        
        Restimate = mean(r,2);
        Rmatrix = toeplitz(Restimate.');    
        [eVecd, eVald] = eig(Rmatrix);
        eVald = diag(eVald);
        [~,sortindexd] = sort(eVald,'descend');
        eVecsortedd = eVecd(:,sortindexd);
        noiseBasisd = eVecsortedd(:,length(us)+1:end);
        Rndirect = noiseBasisd*noiseBasisd';
%%%%%%%%%%%%%%%%%%%NEW ONE
        Restimate = mean(rnew,2);
        Rmatrix = toeplitz(Restimate.');    
        [eVecd, eVald] = eig(Rmatrix);
        eVald = diag(eVald);
        [~,sortindexd] = sort(eVald,'descend');
        eVecsortedd = eVecd(:,sortindexd);
        noiseBasisd = eVecsortedd(:,length(us)+1:end);
        Rndirectnew = noiseBasisd*noiseBasisd';

        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%Algorithm 2 and 3: PRODUCT/MIN MUSIC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        count = 1;
        u = -1:0.001:1;
        Pdirect = zeros(size(u));%%%%direct MUSIC
        Pdirectnew = zeros(size(u));%%%%direct MUSIC
        for idx = u  
            vdirect = exp(-1i*pi * idx*(0:max(lags)).');
            vdirectnew = exp(1i*pi * idx*(0:10).');
            %%%%The three steering vectors above are equal right now, we
            %%%%might change their lengths later
            Pdirect(count) = 1/(vdirect'*Rndirect*vdirect);  
            Pdirectnew(count) = 1/(vdirectnew'*Rndirectnew*vdirectnew);  
            count = count + 1;        
        end
        Pdirect = 10*log10(abs(Pdirect/max(abs(Pdirect))));
        Pdirectnew = 10*log10(abs(Pdirectnew/max(abs(Pdirectnew))));
        if plotfigure
            close all;
            f = figure;    
            ax = axes('Parent', f, 'FontWeight', 'Bold', 'FontSize', 16,... 
            'Position',[0.267203513909224 0.11 0.496339677891654 0.815]);
            hold all;
            plot(ax, u, Pdirect, 'LineWidth', 3, 'Color', [0.6 0 0.6],'LineStyle','-');
            hold on;
            plot(ax, u, Pdirectnew, 'LineWidth', 3, 'Color', [0 0 1],'LineStyle','-.');
            grid on;
            hold on;
            xlabel('u=cos(\theta)', 'FontSize', 16, 'FontWeight', 'Bold');
            ylabel('Output, dB', 'FontSize', 16, 'FontWeight', 'Bold');
            xlim([-1 1]);
            lowerlimit = -15;
            ylim([lowerlimit 0]);
            %%%%The following loop marks the actual source locations by creating a
            %%%%line corresponding to each source location
            for idx = 1:length(us)
                hold on;
                plot([us(idx) us(idx)],[lowerlimit 0],'r:','LineWidth',2);
            end   
            legend('Direct','Partial Direct','Actual u_1','Actual u_2');
            hold on;
            set(gcf,'WindowState','maximized');        
        end
end

function x = ifourierTrans(X,nlower,nhigher,varargin)
%%%%X is the spectrum. nlower is the smallest lag. nhigher is the largest
%%%%lag
    deltau = 0.001;
    if nargin > 3
        deltau = varargin{1};
    end
    u = -1:deltau:1;
    temp1 = length(nlower:nhigher);
    temp2 = size(X,2);
    x = zeros(temp1,temp2);
    count = 1;
    for n = nlower:nhigher
        basis = exp(1i*pi*u*n);
        x(count,:) = 0.5*deltau*basis*X;
        count = count + 1;
    end
end