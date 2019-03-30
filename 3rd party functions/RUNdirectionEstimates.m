function RUNdirectionEstimates
%%%%%For each SNR and number of snapshots, this program calls
%%%%%"directionEstimatesVersion2" 1000 times and saves the MSEs
    close all;clear;clc;
    tic;
    %%%%Design parameters are
    M = 10; N = 12; U1 = 2; U2 = 3;
    snr = 10;
    sdelta = 5;
    snapshots = 20:sdelta:200;
    reps = 1000;
    %%%%We will first run 1000 trials for a fixed snr but a range of snapshots
    for jdx = 1:length(snapshots)
        disp(snapshots(jdx));%%%%Display: Just to get an idea of where we are in the long simulation
        p = zeros(1,reps);    m = zeros(1,reps);    d = zeros(1,reps);    f = zeros(1,reps);
        for ndx = 1:reps
            [p(ndx), m(ndx), d(ndx), f(ndx)] = directionEstimatesVersion2(M,N,U1,U2,snr,snapshots(jdx));
        end
        %%%%the following two lines saves MSEs from all 1000 trials in a
        %%%%file
        b = {'dataForSNRp10SS'};
        save([b{1} num2str(snapshots(jdx))],'p','m','d','f','snr','snapshots');
        disp([mean(p) mean(m) mean(d) mean(f)]);
    end
    toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    close all;clear;clc;
    tic;
    %%%%Design parameters are
    M = 10; N = 12; U1 = 2; U2 = 3;
    snr = -15:15;
    snapshots = 50;
    reps = 1000;
    %%%%Now we will run 1000 trials for a fixed number of snapshots, but
    %%%%different snrs
    for jdx = 1:length(snr)
        disp(snr(jdx));%%%%Display: Just to get an idea of where we are in the long simulation
        p = zeros(1,reps);    m = zeros(1,reps);    d = zeros(1,reps);    f = zeros(1,reps);
        for ndx = 1:reps
            [p(ndx), m(ndx), d(ndx), f(ndx)] = directionEstimatesVersion2(M,N,U1,U2,snr(jdx),snapshots);
        end
        %%%%Save all 1000 MSEs for each algorithm in a file
        b = {'dataForSS50SNRp15Bias'};
        save([b{1} num2str(snr(jdx)+15)],'p','m','d','f','snr','snapshots');
        disp([mean(p) mean(m) mean(d) mean(f)]);
    end
    toc;
end
 