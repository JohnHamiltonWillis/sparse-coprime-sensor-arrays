
snr = [-20,-15,-10,-5,0];
snapshots = 1000;
reps = 10e3;
fileloc = pwd;
filename = 'Direct_vs_PartialMSE_10e3reps.mat';

snr_snapshots_analysis_par2(snr,snapshots,reps,[fileloc,'\',filename]);

for cnt = 1:2
    hold on;
    grid on;
    plot(snr,mse(:,cnt),'-o');
end
xlabel('SNR dB', 'FontSize', 16, 'FontWeight', 'Bold');
ylabel('MSE', 'FontSize', 16, 'FontWeight', 'Bold');
legend('Direct','Partial Direct');
            