function [Plot,mF,w] = MUSICcomp(R,us)
w = linspace(-pi,pi,1000);
P = zeros(1000);
[V,L] = eig(R,'vector');
[sortedL, indices] = sort(L,'descend');
sortedV = V(:,indices);
B = sortedV(:,1+length(us):size(R,2));
for k = 1:1000
    U = exp(1i*w(k)*(0:(size(R,2)-1))');
    D = ctranspose(U)*B;
    temp = D*(ctranspose(D));
    P(k) = 1/temp;
end
P = P/max(abs(P));
mF = 20*log10(abs(P));
w = linspace(-1,1,length(mF));
%%%%Plot
Plot = figure('windowstyle','normal');
axes1 = axes('Parent',Plot,'Position',[0.13 0.11 0.775 0.331]);
hold(axes1,'on');
plot(w,mF,'r','Parent',axes1,'LineWidth',2);
for idx = 1:length(us)
    plot([us(idx) us(idx)],[-60 0],'k:','Parent',axes1,'LineWidth',2);
end
ylabel('Power dB','FontWeight','bold');
xlabel('cosd(\theta)','FontWeight','bold');
title('Spatial Spectral Estimation(CSA)','FontWeight','bold');
xlim(axes1,[-1 1]);
ylim(axes1,[-20 0]);
box(axes1,'on');
grid(axes1,'on');
set(axes1,'FontSize',16,'FontWeight','bold');

legen1 = legend('Estimate','actual');
end