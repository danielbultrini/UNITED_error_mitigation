function rho = init_rho2(qnbs_ex,qnbs_noisy, f)
    global  EC 

    % rho is a (noisy) state sfter initialization or active qubit reset
    % qnbs_ex are noiseless qubits
    % qnbs_n  are noisy 
    % qnbs_ex and qnbs_n are row vectors  
    % f is the noise rescaling factor
    Q = max([qnbs_ex qnbs_noisy]);
    
    trapped_ion_errormodel1_plus_idle(0)
    rex = reshape(EC('ngs'),2,2);
    trapped_ion_errormodel1_plus_idle(f)
    rn = reshape(EC('ngs'),2,2);
    

    if ismember(1,qnbs_ex) 
        r{1} = rex;
    else
        r{1} = rn;
    end
    ii{1} = [-1 -Q-1];
    
    
    for n = 2:Q
        if ismember(n,qnbs_ex) 
            r{n} = rex;
        else
            r{n} = rn;
        end
         ii{n} = [-n -Q-n];
    end
    
    rho = ncon(r,ii);

  



end
