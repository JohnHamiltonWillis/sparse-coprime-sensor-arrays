function DataPreprocessingTest(x)
SampleSize = length(x(1,:));
ApertureEnd = 63;%%%%The array starts at 0 and ends at 63

%Fourier transform trimming
F = zeros(ApertureEnd+1,SampleSize);
S = zeros(ApertureEnd+1,SampleSize);
Strim = zeros(ApertureEnd+1,SampleSize);
xtrim = zeros(ApertureEnd+1,SampleSize);
for idx = 1:ApertureEnd+1
    F(idx,:) = fft(x(idx,:),SampleSize);
    S(idx,:) = fftshift(F(idx,:));
    Strim(idx,:) = [zeros(1,SampleSize/2) S(idx,SampleSize/2+1:SampleSize)];
    xtrim(idx,:) = ifft(Strim(idx,:));
end
w = linspace(-1, 1, length(S));
A = fftshift(fft(x(1,:)));
A = A/max(abs(A));
subplot(2,1,1)
plot(w,20*log10(abs(A)))
xlim([-0.25, 0.25]);
ylim([-30, 0]);
xlabel('x\pi, frequency (rad/s)')
ylabel('S_x(exp(j\omega)) (dB)')
grid on
B = fft(xtrim(1,:));
B = B/max(abs(B));
subplot(2,1,2)
plot(w,20*log10(abs(B)))
xlim([-0.25, 0.25]);
ylim([-30, 0]);
xlabel('x\pi, frequency (rad/s)')
ylabel('S_x(exp(j\omega)) (dB)')
grid on
end

