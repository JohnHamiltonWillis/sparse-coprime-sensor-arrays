U1_list = 2:7;
U2_list = 3:8;
M_list = floor(63./U1_list)+1; % 63 is our array aperture
N_list = floor(63./U2_list)+1;
snr = -10;
snapshots = 100;
reps = 10e3;
fileloc = pwd;
filename = 'valid_config_tech_array_10e3reps.mat';

snr_snapshots_analysis_par(M_list,N_list,U1_list,U2_list,snr,snapshots,reps,[fileloc,'\',filename]);