% function to run data collection for a given number of qubits,
% depth, seed for the RQC, number of non-Clifford gates, 
% training set size and max number of copies.

% the runs the noisy and exact simulations for the circuit of
% interest and the near-Clifford circuits and saves the data in
% .mat files containing matlab cell array structures.
function [] = VD_CDR_data_exe(Q, p, seed, N, Nts, max_copies,...
    nlsp, restart,noise_scaling,save_density_matrices,folder)
    % input arguments must be strings to run on the cluster,
    % convert from strings to numbers here:
    Q = str2num(Q);
    p = str2num(p);
    seed = str2num(seed);
    N = str2num(N);
    Nts = str2num(Nts);
    max_copies = str2num(max_copies);
    nlsp = str2num(nlsp); % noise level state prep.
    restart = str2num(restart); % if code needs to restart from one 
    %training set number, yes = 1, no = 0.
    % generate the RQC:
    
    % [circ,pos,ang] = RQC(p,Q,seed);
    [circ,pos,ang,qbsm] = maxcut_circ(seed,Q);
    % increase the noise by idenity insertion:
    circ_inc_noise = increase_noise_level(circ, nlsp, seed);
    % generate the near-Clifford training circuits:
    in = load(sprintf('tcirc_Q%dseed%dN%dNts%d.mat',Q,seed,N,Nts));
    circnC = in.algproj;
    for i = 1:length(circnC)
    [ctr(i),ctr2(i)] = compare(circ,circnC{i});
    end
    disp([min(ctr) max(ctr) min(ctr2) max(ctr2)]) 
    % define scale factor to reduce noise:
    errtot = total_error(circ,Q);
    if noise_scaling == '1'
        
         f = 1./errtot;
        
        filename1 = sprintf(append(folder,'coi_data_Q%dp%dMC%dnlsp%dseed%d_wns.mat'),Q,p,...
                                        max_copies, nlsp, seed);
        filename3 = sprintf(append(folder,'rhon_Q%dp%dMC%dnlsp%dseed%d_wns.mat'),Q,p,...
                                        max_copies, nlsp, seed);
        filename4 = sprintf(append(folder,'rhoEx_Q%dp%dMC%dnlsp%dseed%d_wns.mat'),Q,p,...
                                        max_copies, nlsp, seed);
    else
        filename1 = sprintf(append(folder,'coi_data_Q%dp%dMC%dnlsp%dseed%d.mat'),Q,p,...
                                        max_copies, nlsp, seed);
        filename3 = sprintf(append(folder,'rhon_Q%dp%dMC%dnlsp%dseed%d.mat'),Q,p,...
                                        max_copies, nlsp, seed);
        filename4 = sprintf(append(folder,'rhoEx_Q%dp%dMC%dnlsp%dseed%d.mat'),Q,p,...
                                        max_copies, nlsp, seed);
        f = 1;
    end
    % now for each training circuit evalutate the circuit with 
    % a different number of copies:
    % first check if this is the first time the code is run, 
    % if so, calculate the circuit of interest data first:
    if restart == 0
        obsex_array_train = zeros(Nts,1);
        obsn_array_train = zeros( Nts,1);
        prob0_array_train = zeros( Nts, max_copies-1);
        prob0p_array_train = zeros( Nts, max_copies-1);
        % calculating the exact result here is a little
        % inefficient as really the result will be the same
        % for each noise level (previously was not a bottle neck).
        [obsex_coi,rhoexact] = obs_exact(circ,Q,qbsm);
        obsex_coi
        % first calculate the noisy, exact and mitigated results:
        [obsn_coi, prob0_array_coi,...
            prob0p_array_coi,rhonoisy] = obs(circ_inc_noise, Q, max_copies,f,qbsm); 
        % organise data into a cell array:
        coi_data = {obsex_coi, obsn_coi, prob0_array_coi,...
            prob0p_array_coi, f*errtot};
        
        % save the data with aprticulat file name:

        save(filename1, 'coi_data');
        if save_density_matrices == '1'
            save(filename3, 'rhonoisy');
            save(filename4, 'rhoexact');
        end
        c_ = 1;
    % if this code has been restarted then find the restart point
    % and run from there:
    elseif restart == 1
        if noise_scaling == '1'
        filename2 = sprintf(...
            append(folder,'train_data_Q%dp%dN%dNts%dMC%dnlsp%dseed%d_wns.mat'),...
            Q,p,N,Nts,max_copies, nlsp, seed );    
        else 
        filename2 = sprintf(...
            append(folder,'train_data_Q%dp%dN%dNts%dMC%dnlsp%dseed%d.mat'),...
            Q,p,N,Nts,max_copies, nlsp, seed );  
        end
        train_data = load(filename2, 'train_data');
        [obsex_array_train,...
        obsn_array_train,...
        prob0_array_train,...
        prob0p_array_train,...
        c_] = train_data.train_data{1,1:5};   
    end
    % run the near-Clifford circuit data.
    % loop through each circuit in the trianing set:
    for c = c_:length(circnC)
        % make tracker file, (can be useful on big runs).
        filename = sprintf(append(folder,'tracker_seed_%d_nlsp_%d.txt'),...
            seed, nlsp);
        fid = fopen(filename,'w');
        string = sprintf('seed %d, circ %d',seed,c);
        fprintf(fid, string);
        fclose(fid);
        % select training circuit:
        tcirc = circnC{c};
        % increase the noise:
        tcirc_inc_noise = increase_noise_level(tcirc, nlsp, seed);
        % calculate the useful information:
            [obsn_train, prob0_train,...
        prob0p_train]= obs(tcirc_inc_noise, Q, max_copies,f, qbsm);
        obsex = obs_exact(tcirc,Q, qbsm)
        % store the data:
        obsn_array_train( c) = obsn_train;
        obsex_array_train( c) = obsex;
        prob0_array_train(c,:) = prob0_train;
        prob0p_array_train(c,:) = prob0p_train;
        train_data = {obsex_array_train, obsn_array_train,...
            prob0_array_train, prob0p_array_train, c, f*errtot};
        % save the data under a file:
        if noise_scaling == '1'
        filename2 = sprintf(...
            append(folder,'train_data_Q%dp%dN%dNts%dMC%dnlsp%dseed%d_wns.mat'),...
            Q,p,N,Nts,max_copies, nlsp, seed );    
        else 
        filename2 = sprintf(...
            append(folder,'train_data_Q%dp%dN%dNts%dMC%dnlsp%dseed%d.mat'),...
            Q,p,N,Nts,max_copies, nlsp, seed );  
        end
        save(filename2, 'train_data');
    end
end

function [ctrnC,ctragr] = compare(alg,algproj)
eps = 1.e-10;
ctrnC = 0;
ctragr = 0;
for i = 1:size(algproj,1)
   if mod(algproj{i,1},pi/4) < eps ||  mod(algproj{i,1},pi/4) > pi/4-eps
   else 
    ctrnC = ctrnC + 1;
   if abs(algproj{i,1} - alg{i,1}) < eps && strcmp(alg{i,3},'Z')
     ctragr = ctragr+1;
   else
     disp([algproj{i,1} alg{i,1}]);
   end 
   end
end
end 
