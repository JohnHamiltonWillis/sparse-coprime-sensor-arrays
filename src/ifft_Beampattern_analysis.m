%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Senior Design - Min and Product Processing Analysis
% Louisiana Tech University
% Pablo Johnson, Daniel Sartori, Tyler Trosclair, John Willis
% Sponsored by Dr. Kaushallya Adhikari
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 2;
N = 3;
U1 = N;
U2 = M;
add1 = U1*1;
add2 = U2*1;

[Bmin, Bprod, u] = ProductMinBeampattern(M, N, U1, U2, add1, add2);
close gcf;
Array = generateSpacings(M,N,1);

x = ifourierTrans(Bmin',Array.lags(1,1),Array.lags(1,end));
R = toeplitz(x);

[V,L] = eig(R,'vector');
[sortedL, indices] = sort(L,'descend');
sortedV = V(:,indices);
w = linspace(-pi,pi,1000);
P = zeros(size(w));
for k = 1:1000
    u = exp(1j*w(k)*(0:4))';
    D = u'*sortedV;
    temp = D*D';
    P(k) = 1/temp;
end
P = P / max(abs(P));
figure;
plot(w/pi,10*log10(abs(P)));
ylim([-30 0]);
xlabel('x pi, rad/s');
ylabel('PSD, db');
grid on;

