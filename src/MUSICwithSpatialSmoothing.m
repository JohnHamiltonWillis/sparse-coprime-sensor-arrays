%%%%Me is the number of sensors in Subarray 1, Ne is the number of
%%%%sensors in Subarray 2, U1 is the undersampling factor of Subarray
%%%%1, U2 is the undersampling factor of Subarray 2, SampleSize is the
%%%%number of snapshots, SNRdB is the SNR in dB. The variable
%%%%plotfigure determines whether the program plots the graphs or not.
%%%%This program is different from the previous version
%%%%"directionEstimates in that it doesn't have "plotfigure" argument
%%%%For each trial, this program keeps passing through a loop until a
%%%%satisfactory data set is obtained
clear all;
close all;
N=12;
M=12;
U1 = 2;
U2 = 3;
num = 3; %%number of sources
SNRdB = 5;
SampleSize = 1000;
us = cosd(randi(181,[1 num])-1);%%%Directions are uniformly distributed from 0 to 180 degrees
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
ApertureEnd = 63;%%%%The array starts at 0 and ends at 63

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

%%%%%%%%%%%%%%%%%%%%ALGORITHM1: Direct MUSIC%%%%%%%%%%%%%%%%%%
%%%%%%Algorithm1 is direct MUSIC. The covariance estimate at each
%%%%%%lag is obtained by taking the average of all approprirate
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
%%%%covariance estimates
r = zeros(length(coarray),length(coarray));
for kdx = 1:SampleSize
    dataset = xtotal(:,kdx);
    %%%%Estimate covariance for a single time sample
    tempR = dataset.*conj(transpose(dataset));
    r = r + tempR;
end
Restimate = r/SampleSize;
%%%%%%% Spatial Smoothing Step:
%%%%%%%Vectorize Restimate and average repeated lags. We may think of z1 as data
%%%%%%%received by a new ULA with sensors from 0 to U1*U2.
z1 = zeros(1,length(coarray));
lags = zeros(1,length(coarray));
%%%%Search for kth lag and store in z1(k). Searches through all columns in
%%%%row 1, then searches through all columns after column 1 in row 2, then
%%%%searches through all columns after column 2 in row 3, ...etc.
for i = 1:length(coarray)
    for j = i:length(coarray)
          k = j-i+1;
            if z1(k)==0
                if abs(Restimate(i,j)) > 0
                    z1(k) = z1(k)+Restimate(i,j); %%sum of each R(k)
                    lags(k) = 1+lags(k); %%number of lags
                end
            end
       
    end
end
%%%%Dividing by lags to take the average
z1 = z1./lags;
z1 = z1(1:U1*U2+1);
%%%%Create matrix from z1 to apply spatial smoothing.
z1matrix = toeplitz(z1);
%%%%% This step below is taking a covariance estimation of z1matrix.
%%%%% z1 should be continuous
z1indicator = ones(1,U2*U1+1); 
%%%%% number of lags produced by z1
coarray2 = conv(z1indicator.',fliplr(z1indicator.'));
coarray2(1:U1*U2) = [];
%%%%covariance estimates
r = zeros(length(coarray2),length(coarray2));
for kdx = 1:U1*U2+1
    dataset = z1matrix(:,kdx);
    %%%%Estimate covariance for a single sensor
    tempR = dataset.*conj(transpose(dataset));
    r = r + tempR;
end
Restimatess = r/(U1*U2+1);
[Plot,mF,w] = MUSICcomp(Restimatess,us);
