M = 10; N = 12; U1 = 2; U2 = 3;

snr = -20:0.1:20;
snapshots = 20:5:200;
mse = zeros(length(snr),length(snapshots),4); %need 4 for p,m,d,f mse's
reps = 1000; % 1000 reps for each snr+snapshot var
params_mse = zeros(reps,4); 
flags_k = zeros(reps,4);
flags = mse; % needs same preallocation

for i_snr = snr
    for j_snapshot = snapshots
        for k = 1:reps
%             [p(k),m(k),d(k),f(k),flag_cnt(k)] = ...
            [params_mse(k,:),flags_k(k,:)] = ...
                directionEstimatesVersion2(M,N,U1,U2,i_snr,j_snapshot);           
        end
        mse(i_snr,j_snapshot,:) = average(params_mse);
        flags = average(flags_k);
    end
end
save([pwd, '\snr_snapshot_mse.mat']);
