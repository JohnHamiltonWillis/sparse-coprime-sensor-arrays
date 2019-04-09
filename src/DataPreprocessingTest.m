clear
SNRdB = -10; SampleSize = 1e3;
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
v = zeros(ApertureEnd+1,numSources);
%%%%%The vector will have the data for all sensors
x = zeros(ApertureEnd+1,SampleSize);
for idx = 1:numSources
    v(:,idx) = exp(1i*kx(idx)*(0:ApertureEnd)*d).';
    x = x + v(:,idx)*s(idx,:);
end
%%%Add proper Gaussian noise samples
x = x + sqrt(varn/2)*randn(ApertureEnd+1,SampleSize) + 1i*sqrt(varn/2)*randn(ApertureEnd+1,SampleSize);

%Fourier transform trimming
F = zeros(ApertureEnd+1,SampleSize);
S = zeros(ApertureEnd+1,SampleSize);
Strim = zeros(ApertureEnd+1,SampleSize);
xtrim = zeros(ApertureEnd+1,SampleSize);
for idx = 1:ApertureEnd+1
    F(idx,:) = fft(x(idx,:));
    S(idx,:) = fftshift(F(idx,:));
    Strim(idx,:) = [zeros(1,SampleSize/2) S(idx,SampleSize/2+1:SampleSize)];
    xtrim(idx,:) = ifft(Strim(idx,:));
end

