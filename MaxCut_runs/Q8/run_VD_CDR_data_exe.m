% function to run data collection for a given number of qubits,
% depth, seed for the RQC, number of non-Clifford gates, 
% training set size and max number of copies.

% the runs the noisy and exact simulations for the circuit of
% interest and the near-Clifford circuits and saves the data in
% .mat files containing matlab cell array structures.
function [] = run_VD_CDR_data_exe(Qt, pt, seed, N, Nts, max_copies,...
    nlsp, restart,noise_scaling,save_density_matrices,folder)
clear global
Qt = str2num(Qt);
pt = str2num(pt);

for ctr = 1:length(Qt)
    Q = Qt(ctr);
    p = pt(ctr);
    VD_CDR_data_exe(num2str(Q), num2str(p), seed, N, Nts, max_copies,...
    nlsp, restart,noise_scaling,save_density_matrices,folder)
end

   
end
