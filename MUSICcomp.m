%%%%%%%%%%Apply MUSIC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mF,w] = MUSICcomp(W,num)
W = conj(W);
R = toeplitz(W);
w = linspace(-pi,pi,1000);
P = zeros(1000);
[V,L] = eig(R,'vector');
[sortedL, indices] = sort(L,'descend');
sortedV = V(:,indices);
B = sortedV(:,1+num:size(W,2));
for k = 1:1000
    U = exp(1i*w(k)*(0:(size(W,2)-1))');
    D = ctranspose(U)*B;
    temp = D*(ctranspose(D));
    P(k) = 1/temp;
end
P = P/max(abs(P));
mF = 20*log10(abs(P));
w = linspace(-1,1,length(mF));