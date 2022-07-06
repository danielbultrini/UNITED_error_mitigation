% This is a script which shows how to execute the code which collects 
% the data neccesary to do VD, CDR and vnCDR and their various
% combinations. The data collection code could be made a little more
% efficient, but I think its okay for now. Feel free to modify. 

% lets just run some data with a 4 qubit system, 3 noise levels up to 6 
% copies and and 20 training circuits with depth 5.
% For now we will just use one RQC seed. We are just measuring a sigma z
% on the edge qubit of our system

% add some paths
path(pathdef)
addpath('RQC_code')
addpath('Data')
addpath('VD_code')
addpath('VD_code_perfect_eff')
addpath('ion_trap_simulator')
addpath('execution_code')
addpath('Data')

%define variables:
Q = '6'; % number of qubits.
p = '5'; % depth of RQC.
seed = '0'; % seed used to generate RQC.
N = '10'; % number of non-Clifford gates.
Nts = '100'; % number of training circuits. 
max_copies = '6'; % max number of copies in VD.
restart = '0';

% loop through the noise levels:
tic
for Q = [4]
    for seed = [0]
       disp(seed)
       for nlsp = [1]
           seed = num2str(seed);
           Q = num2str(Q);
           nlsp_str = num2str(nlsp);
           disp(nlsp_str)
           VD_CDR_data_exe(Q, p, seed, N, Nts, max_copies, nlsp_str, restart,'0','1','')
       end
    end
end
toc

% this code runs in around 2 mins on my computer and generates sample data 
% which we can use to do some VD/CDR comparison. 
% we also want to explore shot noise and how this effects the results. 
% In the data we obtain in VD we measure p0 and p0p, so these are the
% quantaties which will be affected by shot noise. 
% the add_shot_noise function reads the data in the files specified by the
% arguments and then adds simulated shot noise to each noisy observable.
% (the exact values are left untouched).
% shots = 10^5;
% tic
% for Q=[4]
%     for seed = [0]
%         for nlsp = [1]
%             seed = num2str(seed);
%             Q = num2str(Q);
%             add_shot_noise_to_data(shots, str2num(Q), str2num(p), str2num(seed),...
%                  str2num(N), str2num(Nts), str2num(max_copies), nlsp,'0','1','');
%         end
%     end 
% end
% toc

% in the files with 'coi' in the name this means circuit of interest.
% the data is sorted into columns: obs_exact, obs_noisy, p0_array, 
% p0p_array, where obsm = (2*prob0p-1)/(2*prob0-1).

% the files with 'train_data' have the same structure but with with a
% highe dimension becuase they contain more cirucits.

