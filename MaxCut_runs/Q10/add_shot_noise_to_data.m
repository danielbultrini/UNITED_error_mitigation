function [] = add_shot_noise_to_data(shots, Q, p, seed, N, Nts, max_copies, nlsp,noise_rescaling,folder,tag)
    % first we load the data:
    if noise_rescaling == '1'
        filename1 = sprintf(append(folder,'coi_data_Q%dp%dMC%dnlsp%dseed%d_wns.mat'),Q,p,...
                                            max_copies, nlsp, seed);
        coi_data_2 = load(filename1, 'coi_data');
        filename2 = sprintf(...
                append(folder,'train_data_Q%dp%dN%dNts%dMC%dnlsp%dseed%d_wns.mat'),...
                Q,p,N,Nts,max_copies, nlsp, seed );                       
        train_data_2 = load(filename2, 'train_data');
    else
        filename1 = sprintf(append(folder,'coi_data_Q%dp%dMC%dnlsp%dseed%d%s.mat'),Q,p,...
                                            max_copies, nlsp, seed,tag);
        coi_data_2 = load(filename1, 'coi_data');
        filename2 = sprintf(...
                append(folder,'train_data_Q%dp%dN%dNts%dMC%dnlsp%dseed%d%s.mat'),...
                Q,p,N,Nts,max_copies, nlsp, seed,tag );                       
        train_data_2 = load(filename2, 'train_data');
    end
    coi_data_2 = coi_data_2.coi_data;
    train_data_2 = train_data_2.train_data;
    % then we will loop through and add shot noise.
    % first do coi_data:
    coi_obsn_wsn = add_shot_noise(coi_data_2{1,2}, shots);
    prob0_array_coi = arrayfun(@(x) add_shot_noise(x, shots), coi_data_2{1,3});
    prob0p_array_coi = arrayfun(@(x) add_shot_noise(x, shots), coi_data_2{1,4});
    if noise_rescaling == '1'
        filename1_ = sprintf(append(folder,'coi_data_Q%dp%dMC%dnlsp%dseed%dshots%d_wns.mat'),Q,p,...
                                        max_copies, nlsp, seed, shots);
        filename2_ = sprintf(...
            append(folder,'train_data_Q%dp%dN%dNts%dMC%dnlsp%dseed%dshots%d_wns.mat'),...
            Q,p,N,Nts,max_copies, nlsp, seed, shots );
    else 
        filename1_ = sprintf(append(folder,'coi_data_Q%dp%dMC%dnlsp%dseed%dshots%d%s.mat'),Q,p,...
                                        max_copies, nlsp, seed, shots,tag);
        filename2_ = sprintf(...
            append(folder,'train_data_Q%dp%dN%dNts%dMC%dnlsp%dseed%dshots%d%s.mat'),...
            Q,p,N,Nts,max_copies, nlsp, seed, shots,tag );  
    end
    coi_data = {coi_data_2{1,1}, coi_obsn_wsn, prob0_array_coi,...
            prob0p_array_coi};
    save(filename1_, 'coi_data')
    % now do the same thing for the training data:
    train_obsn_wsn = arrayfun(@(x) add_shot_noise(x, shots), train_data_2{:,2});
    prob0_array_train = arrayfun(@(x) add_shot_noise(x, shots), train_data_2{:,3});
    prob0p_array_train = arrayfun(@(x) add_shot_noise(x, shots), train_data_2{:,4});
    
    train_data = {train_data_2{:,1}, train_obsn_wsn, prob0_array_train,...
            prob0p_array_train};
    save(filename2_, 'train_data')

   