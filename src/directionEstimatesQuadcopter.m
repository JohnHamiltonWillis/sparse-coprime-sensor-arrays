function directionEstimatesQuadcopter(M, N, U1, U2,SampleRange)

%WARNING: This script opens all plots at once after completion. It may bog
%down your graphics or crash MATLAB if there are more than 30-50 figures.
    uiopen('*.mat');
    sz = size(totalData);
    sz = sz(1);%The number of samples taken
    
    StartingSample=1;%starting sample desired
    EndingSample=50000;%ending Sample desired
    

    if EndingSample < SampleRange
        EndingSample = StartingSample;
    end
    if EndingSample >= (sz-SampleRange)
        EndingSample = sz - SampleRange;
    end
     %Valid starting values for SampleRange are 1 or SampleRange
    if StartingSample < SampleRange
        StartingSample = SampleRange;
    end

        
    Start = StartingSample/SampleRange;
    End = EndingSample/SampleRange;
    
    %%%%M is the number of sensors in Subarray 1, N is the number of
    %%%%sensors in Subarray 2, U1 is the undersampling factor of Subarray
    %%%%1, U2 is the undersampling factor of Subarray 2, SampleRange creates the
    %%%%range of samples to be used per plot
    
    
    us = cosd([118 62]); %is set to edges                  
    %These parameters are used in the steering vector v
    lambda = 99999999999999999999999;% unused    
    d = lambda/2;    kx = 2*pi/lambda * us;
    
    ApertureEnd = 63;%%%%The array starts at 0 and ends at 63
    
    totalData_r = Rearrange_Data(totalData);               
    for qdx = Start:End
        
        RealSampleSize = (qdx*SampleRange):(qdx*SampleRange+SampleRange-1);%Sets the range of samples to be used
        
        x = totalData_r(RealSampleSize,:); %The feild data within the sample range is called to x

        %Lambda must be defined when using vTotal in min and prod processing
        %aglorithm, and therefore d. These are the only algorithms within
        %directionEstimates that rely on these array parameters. We might need
        %a function that can estimate the frequency even to a relative accuracy
        %so that the value for lambda is accurate. This is to make the script
        %more generalized rather than just manually entering the lamda for the
        %frequency we fired at the array.

        x = x.';%dimensions are flipped in our real data

        indexa = (0:U1:(M-1)*U1).';    indexb = (0:U2:(N-1)*U2).';
        indexunion = unique([indexa' indexb']);
        %%%xa will be subarray 1 data. It will have zero values where
        %%%sensors are not present.
        xa = zeros(max(indexa)+1,length(RealSampleSize));
        %%%xb will be subarray 2 data. It will have zero values where
        %%%sensors are not present
        xb = zeros(max(indexb)+1,length(RealSampleSize));%%%This will be subarray 2 data
        xa(indexa+1,:) = x(indexa+1,:);%%%xa takes the data from x, but only where Subarray 1 has sensors
        xb(indexb+1,:) = x(indexb+1,:);%%%xb takes the data from x, but only where Subarray 2 has sensors
        xtotal(indexunion+1,:) = x(indexunion+1,:);%%%%this is the union of Subarray 1 and Subarray 2. This data will be
                                                   %%%%used in direct MUSIC

        %%%%%%%%%%%%%%%%%%%%ALGORITHM1: Direct MUSIC%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%Algorithm1 is direct MUSIC. The covariance estimate at each
        %%%%%%lag is obtained by taking the average of all approximate
        %%%%%%sensor combinations. The convolution operation  comes in very
        %%%%%%handy in evaluating covariance estimates using the right
        %%%%%%sensor pairs as shown below:
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
        r = zeros(length(coarray),SampleRange);%%%%covariance estimates
        for kdx = 1:SampleRange
            dataset = xtotal(:,kdx); % kdx instead of 1 ?
            %%%%The convolution operation can actually be used to find
            %%%%autocorrelation as shown below, for each set of samples
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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%Algorithm 2 and 3: PRODUCT/MIN MUSIC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        deltau = 0.001; %May increase as needed
        u = -1:deltau:1;
        yprod = zeros(size(u));
        ymin = zeros(size(u));
        %%%%%Apply product/min processing first
        for idx = 1:length(u)
            totalv = (exp(1i*2*pi/lambda * u(idx) *(0:(ApertureEnd+1)).'*d)); % what is this?
            wa = zeros(max(indexa)+1,1);        
            wb = zeros(max(indexb)+1,1);        
            %%array weights
            wa(indexa+1,:) = totalv(indexa+1,:)/M;
            wb(indexb+1,:) = totalv(indexb+1,:)/N;
            tempa = wa'*xa;    
            tempb = wb'*xb;    
            ymin(idx) = sum(min(abs([tempa;tempb])))/SampleRange;
            yprod(idx) = sum(tempa.*conj(tempb))/SampleRange;
        end
        yprod = yprod/max(abs(yprod));
        ymin = ymin/max(abs(ymin));
        %%%%%Eigen values and vectors for product first
        Restimate = ifourierTrans(yprod.', 0,max(lags));
        Rmatrix = toeplitz(Restimate.');
        [eVecd, eVald] = eig(Rmatrix);
        eVald = diag(eVald);
        [~,sortindexd] = sort(eVald,'descend');
        eVecsortedd = eVecd(:,sortindexd);
        noiseBasisd = eVecsortedd(:,length(us)+1:end);
        Rnprod = noiseBasisd*noiseBasisd';
        %%%%%Eigen values and vectors for min next
        Restimate = ifourierTrans((ymin.').^2,0,max(lags));%%%%Remember to square y for min
        Rmatrix = toeplitz(Restimate.');
        [eVecd, eVald] = eig(Rmatrix);
        eVald = diag(eVald);
        [~,sortindexd] = sort(eVald,'descend');
        eVecsortedd = eVecd(:,sortindexd);
        noiseBasisd = eVecsortedd(:,length(us)+1:end);
        Rnmin = noiseBasisd*noiseBasisd';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%MUSIC with FULL ULA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %     indexf = 0:ApertureEnd;    
    %     indexf = indexf';    
    %     xf = zeros(max(indexf)+1,SampleRange);   
    %     xf(indexf+1,:) = x(indexf+1,:);   
    %     Rf = xf*(xf')/SampleRange;
    %     [eVec, eVal] = eig(Rf);
    %     eVal = diag(eVal);
    %     [~,sortindex] = sort(eVal,'descend');
    %     eVecsorted = eVec(:,sortindex);
    %     noiseBasis = eVecsorted(:,length(us)+1:end);
    %     Rnf = noiseBasis*noiseBasis';    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        count = 1;
        Pdirect = zeros(size(u));%%%%direct MUSIC
        Pmin = zeros(size(u));%%%%Min MUSIC
        Pprod = zeros(size(u));%%%%Product MUSIC
    %     Pf = zeros(size(u));%%%%MUSIC for full ULA
        for idx = u  
            vprod = exp(-1i*pi * idx*(0:max(lags)).');
            vmin = exp(-1i*pi * idx*(0:max(lags)).');
            vdirect = exp(-1i*pi * idx*(0:max(lags)).');
            %%%%The three steering vectors above are equal right now, we
            %%%%might change their lengths later

    %         vf = exp(1i*pi * idx*(0:ApertureEnd).');
            Pprod(count) = 1/(vprod'*Rnprod*vprod);
            Pmin(count) = 1/(vmin'*Rnmin*vmin);
            Pdirect(count) = 1/(vdirect'*Rndirect*vdirect);    
    %         Pf(count) = 1/(vf'*Rnf*vf);
            count = count + 1;        
        end
        %%%%%Normalize and convert to dB
        Pprod = 10*log10(abs(Pprod/max(abs(Pprod))));
        Pmin = 10*log10(abs(Pmin/max(abs(Pmin))));
        Pdirect = 10*log10(abs(Pdirect/max(abs(Pdirect))));
    %     Pf = 10*log10(abs(Pf/max(abs(Pf))));
        %%%%%%The figure with direction estimates are plotted in every
        %%%%%%trial, but we won't see it because of "close all". But when 
        %%%%%%are debugging with breakpoints, we might want to see the
        %%%%%%figures
        f = figure;    
        ax = axes('Parent', f, 'FontWeight', 'Bold', 'FontSize', 16,... 
        'Position',[0.267203513909224 0.11 0.496339677891654 0.815]);
        hold all;
        plot(ax, u, Pprod, 'LineWidth', 3, 'Color', 'Black','LineStyle', '-');
        hold on;
        plot(ax, u, Pmin, 'LineWidth',3, 'Color', 'b','LineStyle', '--');
        hold on;
        plot(ax, u, Pdirect, 'LineWidth', 3, 'Color', [0.6 0 0.6],'LineStyle','-.');
        hold on;
%         plot(ax, u, Pf, 'LineWidth', 3, 'Color', [0 0.6 0], 'LineStyle', ':');
%         grid on;
%         hold on;
        xlabel('u=cos(\theta)', 'FontSize', 16, 'FontWeight', 'Bold');
        ylabel('Output, dB', 'FontSize', 16, 'FontWeight', 'Bold');
        xlim([-1 1]);
        lowerlimit = -20;
        ylim([lowerlimit 0]);
        %%%The following loop marks the actual source locations by creating a
        %%%line corresponding to each source location
        for idx = 1:length(us)
            hold on;
            plot([us(idx) us(idx)],[lowerlimit 0],'r:','LineWidth',2);
        end   
        legend('Product','Min','Direct');
        hold on;
        set(gcf,'WindowState','maximized');
    end
%%%%%The rest of the program finds the peaks in our estimates and
%%%%%computes the Mean Squared Errors
%     MinPeakHeight = -12;
%    [~,prod_locs] = findpeaks(Pprod,u,'NPeaks',2,'MinPeakHeight',MinPeakHeight);
%    [~,min_locs] = findpeaks(Pmin,u,'NPeaks',2,'MinPeakHeight',MinPeakHeight);
%    [~,direct_locs] = findpeaks(Pdirect,u,'NPeaks',2,'MinPeakHeight',MinPeakHeight);
%    [~,full_locs] = findpeaks(Pf,u,'NPeaks',2,'MinPeakHeight',MinPeakHeight);
     %%%%Compute the MSE. We don't know which peak locations correspond
     %%%%with which directions. So, we will associate _locs(1) with us(1)
     %%%%and compute the total MSE and call it mse1. Then, we will
     %%%%associate _locs(1) with us(2) and compute the total MSE and call
     %%%%it mse2. Then, actual MSE = min(mse1,mse2). Also, the
     %%%%estimate might have only one peaks. To account for that, first
     %%%%check the length of _locs and _locs. If they are not length 2,
     %%%%make them length 2 by repeating the same peak. end

    %product MSE calculation
%     pMSE1 = sum((us-prod_locs).^2)/2;
%     pMSE2 = sum((fliplr(us)-prod_locs).^2)/2;
%     pMSE = min(pMSE1,pMSE2);    
%     %minimum MSE calculation
%     mMSE1 = sum((us-min_locs).^2)/2;
%     mMSE2 = sum((fliplr(us)-min_locs).^2)/2;
%     mMSE = min(mMSE1,mMSE2);
%     %direct MSE calculation
%     dMSE1 = sum((us-direct_locs).^2)/2;
%     dMSE2 = sum((fliplr(us)-direct_locs).^2)/2;
%     dMSE = min(dMSE1,dMSE2);
    %full MSE calculation
%     fMSE1 = sum((us-full_locs).^2)/2;
%     fMSE2 = sum((fliplr(us)-full_locs).^2)/2;
%     fMSE = min(fMSE1,fMSE2);

end

function [rearranged_data, correct_order] = Rearrange_Data(totalData)
% Rearrange_Data Rearrange data from LaTech Senior Design array so that the
% data is linear from left to right
% second output argument is the correct_order vector. 
correct_order = [16	32	15	31	14	30	13	29	12	28	11	27	10	26	9	25	8	24	7	23	6	22	5	21	4	20	3	19	2	18	1	17	33	49	34	50	35	51	36	52	37	53	38	54	39	55	40	56	41	57	42	58	43	59	44	60	45	61	46	62	47	63	48	64];

rearranged_data = totalData(:,correct_order);
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