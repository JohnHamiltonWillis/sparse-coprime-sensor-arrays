%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Senior Design - MUSIC analysis |comparison between autocorrelation|
%                                |estimators                        |
% Louisiana Tech University
% Pablo Johnson, Daniel Sartori, Tyler Trosclair, John Willis
% Sponsored by Dr. Kaushallya Adhikari
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Function for comparing spatial spectral estimation with MUSIC using two
%%different methods for estimating the autocorrelation matrix. These
%%estimations are applied to both a ULA and CSA and their PSDs are plotted
%%next to eachother for comparison.
%Number of sensors 'A' in ULA, spacings of subarrays N and M (M > N)
%num is the number of signals.
%To compare the two autocorrelation estimations set iteration to
%~1000. Set iteration to 1 to just a single iteration
%(THERE IS NO PLOTTING FOR ITERATION > 1).
function [T,figure1] = MUSIC_comparison(A,N,M,num,iterations)
bothsuccessfulnum = 0;
singlesuccessfulnum = 0;
neithersuccessfulnum = 0;
averagesuccessfulnum = 0;
for sr = 1:iterations             %%remove this for loop and uncomment plot sections
    % to plot rather than table of success rates
    %number of signals
    %%%%%%%%%% Create 'num' signals 'x' buried in noise to be measured %%%
    %Some constant parameters
    angle =  10+(160).*rand(num,1);         %%column vector of random signal
    %directions in degrees
    f = 8923.8;                             %%Signal frequency in Hz
    c = 343;                                %%Signal speed in m/s
    deltat = 1/500;                         %%Temporal sampling interval (s)
    SNRdB = -10;                            %%SNR in dB
    vars = 1;                               %%Signal variance
    varn = vars*10^(-SNRdB/10);             %%Noise variance
    times = (0:deltat:1)';
    locations = 0:(A-1);
    [indices,t] = meshgrid(locations,times);
    sumx = 0;
    %Create a sum 'sumx' of 'num' different signals with random angles
    for i = 1:num
        x = exp(1j*(2*pi*f*t-pi.*cosd(angle(i)).*indices));
        sumx = sumx + x;
    end
    clear indices t times;
    totalData = sumx + sqrt(varn/2)*randn(size(x)) + 1i*sqrt(varn/2)*randn(size(x));
    for array = 1:2
        %%%%%%%%%%%%%%%%%%%%%% CSA data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if array == 1
            %%%%%%%%%% Create index that matches CSA sensor locations %%%%%%%
            %Indices that matches each of the two subarrays
            sparseindices1 = (0:N:(A-1)) + 1;
            sparseindices2 = (0:M:(A-1)) + 1;
            %data collected by the two subarrays
            data1 = zeros(size(totalData));
            data2 = zeros(size(totalData));
            data1(:,sparseindices1) = totalData(:,sparseindices1);
            data2(:,sparseindices2) = totalData(:,sparseindices2);
            %Find union between the column data in two subarrays (Equivalent to CSA data)
            data = zeros(size(totalData));          %data in CSA
            for i = 1:A
                if (data1(:,i) == 0)
                    data(:,i) = data2(:,i);
                elseif (data2(:,i) == 0)
                    data(:,i) = data1(:,i);
                else
                    data(:,i) = data1(:,i);
                end
            end
        end
        %%%%%%%%%%% ULA data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if array == 2
            data = totalData;
        end
        
        %%%%%%%%%% Autocorrelation Estimations %%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Estimation using a single product
        W = [];                                 %Autocorrelation Matrix
        for k = 1:A
            D = [];                             %Product Matrix
            
            %Take the product between two sensors k sensors apart from eachother.
            %One sensor is located at i, the other at i+k-1.
            %i.e. for i=1 the first product in the loop is between first column
            %   in data and a column spaced k apart and each iteration shifts
            %   ...each column over one sensor until a nonzero product is produced.
            
            for i = 1:(A-k+1)
                if abs(data(:,i).*conj(data(:,i+k-1))) > 0
                    D = data(:,i).*conj(data(:,i+k-1));
                    break
                else
                    D = 0;
                end
            end
            if D == 0
                break
            end
            W = [W D];
        end
        W = mean(W);
        [mF1,w] = MUSICcomp(W,num);
        %%% Estimation by averaging several products
        W = [];
        for k = 1:A
            D = [];                             %Product Matrix
            
            %Take the product between two sensors k sensors apart from eachother.
            %One sensor is located at i, the other at i+k-1.
            %i.e. for i=1 the first product in the loop is between first column
            %   in data and a column spaced k apart and each iteration shifts
            %   ...each column over one sensor.
            
            for i = 1:(A-k+1)
                %%Remark: Here it seems that the product for lag(0) is computed
                %  when k=1.
                D = [D data(:,i).*conj(data(:,i+k-1))];
            end
            %Remove empty columns(only interested in products between active mics)
            %Rather than removing zeros, keep them and truncate up to the first
            %zero.
            D = mean(D);
            D = unique(D);
            D = D(2:size(D,2));
            %Prevent division by 0 in mean() (this happens when i=A-k+1 ==> D = [])
            if abs(D) > 0
                Ri = mean(D,2);
                W = [W Ri];
            else
                break
            end
            
        end
        [mF2,w] = MUSICcomp(W,num);
        
        [pks1, locs1] = findpeaks(mF1,w,'MinPeakHeight',-25);
        [pks2, locs2] = findpeaks(mF2,w,'MinPeakHeight',-25);
        if iterations == 1
            %%%%%%%%Plot the spatial spectrum (CSA)
            if array == 1
                figure1 = figure('windowstyle','normal');
                axes1 = axes('Parent',figure1,...
                    'Position',[0.13 0.11 0.775 0.331]);
                hold(axes1,'on');
                plot(w,mF1,'r','Parent',axes1,'LineWidth',2);
                plot(w,mF2,'b','Parent',axes1,'LineWidth',2);
                for idx = 1:num
                    plot([cosd(angle(idx)) cosd(angle(idx))],[-60 0],'k:','Parent',axes1,'LineWidth',2);
                end
                ylabel('Power dB','FontWeight','bold');
                xlabel('cosd(\theta)','FontWeight','bold');
                title('Spatial Spectral Estimation(CSA)','FontWeight','bold');
                xlim(axes1,[-1 1]);
                ylim(axes1,[-20 0]);
                box(axes1,'on');
                grid(axes1,'on');
                set(axes1,'FontSize',16,'FontWeight','bold');
                
                legen1 = legend('single','average');
            end
            if array == 2
                %%%%%%%%Plot the spatial spectrum (ULA)
                axes2 = axes('Parent',figure1,...
                    'Position',[0.13 0.61 0.775 0.331]);
                hold(axes2,'on');
                plot(w,mF1,'r','Parent',axes2,'LineWidth',2);
                plot(w,mF2,'b','Parent',axes2,'LineWidth',2);
                for idx = 1:num
                    plot([cosd(angle(idx)) cosd(angle(idx))],[-60 0],'k:','Parent',axes2,'LineWidth',2);
                end
                ylabel('Power dB','FontWeight','bold');
                xlabel('cosd(\theta)','FontWeight','bold');
                title('Spatial Spectral Estimation (ULA)','FontWeight','bold');
                xlim(axes2,[-1 1]);
                ylim(axes2,[-20 0]);
                box(axes2,'on');
                grid(axes2,'on');
                set(axes2,'FontSize',16,'FontWeight','bold');
                legen2 = legend('single','average');
            end
        end
    end
    %%%%%%%%%%%%%%%%% Success probability comparison %%%%%%%%%%%%%%%%%%%%%
    % counts how often each autocorrelation estimation is successful
    % success is measured simply by the number of peaks - if the number
    %   peaks (locs1 and locs2) matches the number of signals (num) then
    %   signal detection is successful
    % Note: Although this is not a sufficient condition for signal
    %   detection, it is the first condition for signal detection and
    %   and should be a nice approximation for now. This should be made
    %   more precise later.
    if size(locs1,2) == num && size(locs2,2) == num
        bothsuccessfulnum = bothsuccessfulnum + 1;
        
    elseif size(locs2,2) == num
        averagesuccessfulnum = averagesuccessfulnum + 1;
    elseif size(locs1,2) == num
        singlesuccessfulnum = singlesuccessfulnum + 1;
    else
        neithersuccessfulnum = neithersuccessfulnum + 1;
    end
    comparison = [bothsuccessfulnum averagesuccessfulnum singlesuccessfulnum neithersuccessfulnum];
    cr = 100*comparison/sr;  %successrate in percent
    
%John Comparison
    % Calculate the PSD of the sum of signals
    window = chebwin(length(sumx(:,1)), 40);
    windowed_data = window.'.*sumx(:,1)';
    F = fft(windowed_data, length(mF1));
    PSDx = fftshift(F)';

    % Calculate the correlation coefficient between the two
    mF1CC = corrcoef(PSDx,mF1);
    mF1CC(1, 2)
    abs(mF1CC(1, 2))
    mF2CC = corrcoef(PSDx,mF2);
    mF2CC(1, 2)
    abs(mF2CC(1, 2))
end
if iterations > 1
    %%Table display success probabilities in percentage
    T = array2table(cr,'VariableNames',{'both', 'average', 'single', 'neither'})
else
    T = [];
end
