%%%%Me is the number of sensors in Subarray 1, Ne is the number of
%%%%sensors in Subarray 2, U1 is the undersampling factor of Subarray
%%%%1, U2 is the undersampling factor of Subarray 2, SampleSize is the
%%%%number of snapshots, SNRdB is the SNR in dB. The variable
%%%%plotfigure determines whether the program plots the graphs or not.
%%%%This program is different from the previous version
%%%%"directionEstimates in that it doesn't have "plotfigure" argument
%%%%For each trial, this program keeps passing through a loop until a
%%%%satisfactory data set is obtained
N=32;
M=21;
U1 = 2;
U2 = 3;
%%%U1<U2
SNRdB = 10;
SampleSize = 100;
us = cosd(randi(181,[1 2])-1);%%%Directions are uniformly distributed from 0 to 180 degrees
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
%%%%The coarray could have holes. Keep the longest hole-free and
%%%%the associated lags.
% firstzeroindex = find(coarray==0,1);%%%finds the index of the first zero elementconv(a,b)
%%%%if the firstzeroindex is not empty, remove the elements from
%%%%coarray and lag starting at the first zero element
% if ~isempty(firstzeroindex)
%     coarray(firstzeroindex:end) = [];
%     lags(firstzeroindex:end) = [];
% end
%%%%We remove the empty elements later in spatial smoothing instead of
%%%%in the steps above.
r = zeros(length(coarray),SampleSize);%%%%covariance estimates
for kdx = 1:SampleSize
    dataset = xtotal(:,kdx);
    %%%%The convolution operation can actually be used to find
    %%%%autocorrelation as shown below, for each set of samples
    tempR = conv(dataset.',fliplr(conj(dataset.')));
%   tempR(1:temp) = [];
%     if ~isempty(firstzeroindex)
%         tempR(firstzeroindex:end) = [];
%     end
    %%%%The step above is removed for spatial smoothing
    r(:,kdx) = (tempR.')./(coarray.');
end
Restimate = mean(r,2);
Rmatrix = toeplitz(Restimate.');
%%%%%%% Spatial Smoothing Step
%%%%Vectorize Rmatrix and remove repeated lags. We may think of z1 as data
%%%%received by a new ULA with sensors from -U1*U2 to U1*U2.
z1 = zeros(1,U1*N);
for i = 1:U1*N 
    for k = 1:U1*N 
        for j = k:U1*N 
            if z1(k)==0
                if abs(Rmatrix(i,j)) > 0
                    z1(k) = Rmatrix(i,j);
                end
            end
        end
    end
end
%%%%The computation above does not guarantee that the set of negative lags
%%%%in z1 will be continuous. The positive lags are continuous so we will
%%%%use the set of positive lags to create the negative lags.
z1 = z1(1:U1*U2+1);
% z2 = conj(fliplr(z1(2:M*N+1)));
% z3 = [z2, z1];
%%%%Create matrix from z1 to apply spatial smoothing.
z1matrix = toeplitz(z1);
%%%%% This step below is taking a covariance estimation of z1matrix.
zss = zeros(U2*U1+1,U2*U1+1);
%%%%% z1 should be continuous
z1indicator = ones(U2*U1+1); 
%%%%% number of lags produced by z1
coarray2 = conv(z1indicator,fliplr(z1indicator));
for kdx = 1:U2*U1+1
    dataset2 = z1matrix(:,kdx);
    tempR2 = conv(dataset2.',fliplr(conj(dataset2.')));
    tempR2(1:U1*U2) = [];
    zss(:,kdx) = (tempR2.')./(coarray2.');
end
Restimatess = mean(zss,2);
Rmatrixss = toeplitz(Restimatess.');
[eVecd, eVald] = eig(Rmatrixss);
eVald = diag(eVald);
[~,sortindexd] = sort(eVald,'descend');
eVecsortedd = eVecd(:,sortindexd);
noiseBasisd = eVecsortedd(:,length(us)+1:end);
Rndirect = noiseBasisd*noiseBasisd';