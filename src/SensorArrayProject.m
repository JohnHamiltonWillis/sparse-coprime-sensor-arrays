function SensorArrayProject(measurement_angle,sensor_layout,varargin)
m = length(varargin);
OptArgs = {1 25 '../'};
OptArgs(1:m) = varargin;
[flag, nBlocksToGrab, filepath] = OptArgs{:};

if flag == 0
    gatherVDAM(measurement_angle,nBlocksToGrab,filepath);
elseif flag == 1
    uiopen('*.mat');
else
    fprintf('Acceptable flags are 0 to gather data or 1 to use available data\n')
end

cases = CSAFinder(sensor_layout);

for i = 1:length(cases)
    directionEstimatesRealData(totalData,cases(i,1),cases(i,2),cases(i,3),cases(i,4),measurement_angle);
end
end

function gatherVDAM(measurement_angle,nBlocksToGrab,filepath)

block_size = 8192; %this will define the number of samples returned per channel
% nBlocksToGrab = 25; %this will define how many consecutive blocks we wish to record. 


% **** first we create a microphone array object

MA = Microphone_Array();

% **** the next step is to set any parameters, in this case we will set the 
%      block transfer size to 8192 

MA.block_size = block_size;

% **** set the gain setting for the internal analog circuitry, 
%  1 =  -15dB
%  2 =  0dB
%  3 =  15dB
MA.gainSetting = 2;

% **** No we can initialize the array, this will attempt to open the driver
% and connect to the microphone array.  If the array is not connected to
% the computer then this will fail. 

MA.init_array();

% Now we must start the internal acquisition of the microphone array data.
% The driver will begin acquiring data into the background circular buffer
% immediately after this call. 

MA.start();

%At this point we can begin grabbing data from the array.  In order to grab
%data from the array we must first creat a buffer in Matlab to store all
%audio data.  This array should be of size n_channels x block_size

data = zeros(MA.n_channels*MA.block_size,1);

% And we will create a large matrix to store all of the audio data acquired

totalData = zeros(MA.block_size*nBlocksToGrab,MA.n_channels);

% cleanupObj = onCleanup(@() closeVDAM(measurement_angle, totalData,filepath,MA));

% Now we will grab <nBlocksToGrab> consecutive frames 

for i = 1:nBlocksToGrab
    fprintf('Grabbing frame %d\n',i);
    data = getDataBlock(MA,data);
    totalData((i-1)*MA.block_size+1:i*MA.block_size,:) = reshape(data,[MA.block_size,MA.n_channels]);
end

% 
% % now we can write these to n_channel individual wav files
% plot(totalData(:,1:8))
% for i = 1:MA.n_channels
%     fprintf('Writing file number %d\n',i);
%     filename = sprintf('%s_%0.2d.wav',wavPrefix,i);
%     audiowrite(filename, totalData(:,i),MA.sample_rate);
% %     wavwrite(totalData(:,i),MA.sample_rate,filename);
% end

%now we need to clean up the workspace by stopping the acquisition and then
%closing the Microphone_Array object


closeVDAM(measurement_angle,totalData,filepath,MA);


end

function closeVDAM(measurement_angle,totalData,filepath,MA)
    fprintf('saving data, closing VDAM...\n');
    filename = [filepath '\' datestr(now,'yyyy-mm-dd_HHMMSS') '_' measurement_angle, '.mat'];
    save(filename, 'totalData');
    MA.stop();
    MA.close();
    delete(MA);
end

function [coprime_arrays] = CSAFinder(sensor_layout)
% sensor_layout is a vector of 1's and 0's where a 1 indicates that a 
% sensor is available at that position. array_length is the total number of
% sensor positions in the array. This function finds all the possible 
% coprime sensor array layouts that will fit in any given sparse linear 
% array. This returns coprime_arrays which gives M, N, U1, and U2 in each 
% row of the array.
coprime_arrays = [];
array_length = length(sensor_layout);
for spacing = 1:(array_length-1)
    % Iterate through all the coprime pair spacings3 and generate pairs
    cpairs = GenerateCoprimePairs(2,array_length, spacing);
    for pair = 1:length(cpairs)
        % Iterate through all the coprime pairs
        for max_sensor1 = 1:floor((array_length-1)/cpairs{pair}(1)+1)
            % Iterate through all the period extensions that fit in the
            % full sensor array for subarray1
            for max_sensor2 = 1:floor((array_length-1)/cpairs{pair}(2)+1)
                % Iterate through all the period extensions that fit in the
                % full sensor array for subarray2
                % Create the coprime layout necessary for a given coprime
                % pair and their subarrays extensions
                coprime_array = CoprimeArray(max_sensor1,max_sensor2,cpairs{pair}(1),cpairs{pair}(2));
                coprime_array = coprime_array.array;
                for shift = 1:(array_length-length(coprime_array))
                    % Shift coprime array from the beginning to the end of
                    % the full array
                    coprime_array = [0 coprime_array];
                    % Check if the proper sensors are available for coprime
                    % array
                    eligible = true;
                    for sensor = 1:length(coprime_array)             
                        if ((coprime_array(sensor)==1) && (sensor_layout(sensor)==0))
                            eligible = false;
                            break
                        end
                    end
                end
                % Store eligible coprime layout
                if eligible == true
                    coprime_arrays = [coprime_arrays; max_sensor1 max_sensor2 cpairs{pair}(1) cpairs{pair}(2)];
                end
            end
        end
    end
end
end

function [coprimes] = GenerateCoprimePairs(min,max,spacing)
% min is the lower bound in the range being searched. Max is the upper
% bound. Spacing sets the desired difference between the two coprimes. 

if spacing < 1
    spacing = 1;
end
if min < 2 % trivially don't want to test 1 since it's subarray is continuous
    min = 2;
    disp('Min too low, setting min to 2');
end
if min > max % Min must be greater than max
    error('Max must be greater than min');
end
%%
coprimes = cell(1, max); % Preallocate space 
for N = min+spacing:max % higher number n goes from min + spacing to max
    M = N-spacing; % lower number m
    factors_N = unique(factor(N)); % unique factors of N
        
    if ~ismember(M,factors_N) % if M is not a factor of N
        factors_M = unique(factor(M)); % find factors of M
    else
        factors_M = M; % else set M as the only factor of M
    end
    common_factors = intersect(factors_N, factors_M); % find intersection of factors of M and N
    if isempty(common_factors) % if there are no common factors
        coprimes{N} = [M,N]; % save pair as coprime
    end
end
coprimes = coprimes(~cellfun('isempty',coprimes)); % remove empty cell space
end

function Subarray = CoprimeArray(M,N,U1,U2)
% Function to generate vector representation of the coprime array given Me,
% Ne, and the undersamping factors. Also returns the subarrays and the
% available lags and number of each lag. 


max_sensor = max((M-1)*U1, (N-1)*U2)+1;
 
 
% Generate 0 1 representation of subarray 1 and 2
Subarray1 = zeros(1,max_sensor);
Subarray1((0:U1:(M-1)*U1)+1) = 1;
% Subarray1 = [1,Subarray1];
Subarray2 = zeros(1,max_sensor);
Subarray2((0:U2:(N-1)*U2)+1) = 1;
% Subarray2 = [1,Subarray2];

Sensor_placement = max(vertcat(Subarray1, Subarray2)); % Combine two subarrays

% Save the generated data to the Subarray struct
Subarray.array = Sensor_placement;
Subarray.sub1 = Subarray1;
Subarray.sub2 = Subarray2;

% Find the available lags given our coprime array
coarray = conv(Sensor_placement,fliplr(Sensor_placement));
Subarray.coarray = coarray;

end

function [pMSE,mMSE,dMSE,fMSE] = directionEstimatesRealData(totalData, M, N, U1, U2, deg)
    %%%%M is the number of sensors in Subarray 1, N is the number of
    %%%%sensors in Subarray 2, U1 is the undersampling factor of Subarray
    %%%%1, U2 is the undersampling factor of Subarray 2, SampleRange creates the
    %%%%range of samples to be used. In this case SampleRange is multiplied by for the
    %%%%min range and 6 for the max range. The variable
    %%%%plotfigure determines whether the program plots the graphs or not.
    plot_fig = 1;
    
    us = cosd([180-deg deg]);
    %The directions are from 0 to 180 but the results we desire are from 90
                        
    %These parameters are used in the steering vector v
    lambda = 340/3000;% meters %sound speed/frequency;    
    d = lambda/2;    kx = 2*pi/lambda * us;
    
    ApertureEnd = 63;%%%%The array starts at 0 and ends at 63

    RealSampleSize = floor(1/5*length(totalData)):length(totalData);%Sets the range of samples to be used
    SampleRange = length(RealSampleSize);
    
    x = totalData(RealSampleSize,:); %The field data within the sample range is called to x
    
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
    indexf = 0:ApertureEnd;    
    indexf = indexf';    
    xf = zeros(max(indexf)+1,SampleRange);   
    xf(indexf+1,:) = x(indexf+1,:);   
    Rf = xf*(xf')/SampleRange;
    [eVec, eVal] = eig(Rf);
    eVal = diag(eVal);
    [~,sortindex] = sort(eVal,'descend');
    eVecsorted = eVec(:,sortindex);
    noiseBasis = eVecsorted(:,length(us)+1:end);
    Rnf = noiseBasis*noiseBasis';    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    count = 1;
    Pdirect = zeros(size(u));%%%%direct MUSIC
    Pmin = zeros(size(u));%%%%Min MUSIC
    Pprod = zeros(size(u));%%%%Product MUSIC
    Pf = zeros(size(u));%%%%MUSIC for full ULA
    for idx = u  
        vprod = exp(-1i*pi * idx*(0:max(lags)).');
        vmin = exp(-1i*pi * idx*(0:max(lags)).');
        vdirect = exp(-1i*pi * idx*(0:max(lags)).');
        %%%%The three steering vectors above are equal right now, we
        %%%%might change their lengths later

        vf = exp(1i*pi * idx*(0:ApertureEnd).');
        Pprod(count) = 1/(vprod'*Rnprod*vprod);
        Pmin(count) = 1/(vmin'*Rnmin*vmin);
        Pdirect(count) = 1/(vdirect'*Rndirect*vdirect);    
        Pf(count) = 1/(vf'*Rnf*vf);
        count = count + 1;        
    end
    %%%%%Normalize and convert to dB
    Pprod = 10*log10(abs(Pprod/max(abs(Pprod))));
    Pmin = 10*log10(abs(Pmin/max(abs(Pmin))));
    Pdirect = 10*log10(abs(Pdirect/max(abs(Pdirect))));
    Pf = 10*log10(abs(Pf/max(abs(Pf))));
    %%%%%%The figure with direction estimates are plotted in every
    %%%%%%trial, but we won't see it because of "close all". But when 
    %%%%%%are debugging with breakpoints, we might want to see the
    %%%%%%figures
    if plot_fig
        close all;
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
        plot(ax, u, Pf, 'LineWidth', 3, 'Color', [0 0.6 0], 'LineStyle', ':');
        grid on;
        hold on;
        xlabel('u=cos(\theta)', 'FontSize', 16, 'FontWeight', 'Bold');
        ylabel('Output, dB', 'FontSize', 16, 'FontWeight', 'Bold');
        xlim([-1 1]);
        lowerlimit = -20;
        ylim([lowerlimit 0]);
        %%%%The following loop marks the actual source locations by creating a
        %%%%line corresponding to each source location
        for idx = 1:length(us)
            hold on;
            plot([us(idx) us(idx)],[lowerlimit 0],'r:','LineWidth',2);
        end   
        legend('Product','Min','Direct','Full ULA','Actual u_1','Actual u_2');
        title(['M = ',num2str(M),' N = ',num2str(N),' U1 = ',num2str(U1),' U2 = ',num2str(U2)])
        hold on;
        set(gcf,'WindowState','maximized');
        savefig(f,['DirectionEstimates','_',num2str(M),'_',num2str(N),'_',num2str(U1),'_',num2str(U2),'_',num2str(deg),'_',datestr(now,'yyyy-mm-dd-HHMMSS'),'.fig'])
    end

%%%%%The rest of the program finds the peaks in our estimates and
%%%%%computes the Mean Squared Errors
    MinPeakHeight = -12;
   [~,prod_locs] = findpeaks(Pprod,u,'NPeaks',2,'MinPeakHeight',MinPeakHeight);
   [~,min_locs] = findpeaks(Pmin,u,'NPeaks',2,'MinPea kHeight',MinPeakHeight);
   [~,direct_locs] = findpeaks(Pdirect,u,'NPeaks',2,'MinPeakHeight',MinPeakHeight);
   [~,full_locs] = findpeaks(Pf,u,'NPeaks',2,'MinPeakHeight',MinPeakHeight);
     %%%%Compute the MSE. We don't know which peak locations correspond
     %%%%with which directions. So, we will associate _locs(1) with us(1)
     %%%%and compute the total MSE and call it mse1. Then, we will
     %%%%associate _locs(1) with us(2) and compute the total MSE and call
     %%%%it mse2. Then, actual MSE = min(mse1,mse2). Also, the
     %%%%estimate might have only one peaks. To account for that, first
     %%%%check the length of _locs and _locs. If they are not length 2,
     %%%%make them length 2 by repeating the same peak. end

    %product MSE calculation
    pMSE1 = sum((us-prod_locs).^2)/2;
    pMSE2 = sum((fliplr(us)-prod_locs).^2)/2;
    pMSE = min(pMSE1,pMSE2);    
    %minimum MSE calculation
    mMSE1 = sum((us-min_locs).^2)/2;
    mMSE2 = sum((fliplr(us)-min_locs).^2)/2;
    mMSE = min(mMSE1,mMSE2);
    %direct MSE calculation
    dMSE1 = sum((us-direct_locs).^2)/2;
    dMSE2 = sum((fliplr(us)-direct_locs).^2)/2;
    dMSE = min(dMSE1,dMSE2);
    %full MSE calculation
    fMSE1 = sum((us-full_locs).^2)/2;
    fMSE2 = sum((fliplr(us)-full_locs).^2)/2;
    fMSE = min(fMSE1,fMSE2);
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