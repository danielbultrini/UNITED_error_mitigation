function rho = evolve(algc,rho,qnbs_ex,qnbs_n, f)
    global  EC
    % rho is a  (nosiy) state
    % algc is a compiled circuit
    % qnbs_ex are noiseless qubits
    % qnbs_n  are noisy 
    % qnbs_ex and qnbs_n are row vectors  
    % tehere's no checkup of consistancy in between total number of qubits
    % in algc and qnbs_ex and qnbs_n
    
    Q = max([qnbs_ex qnbs_n]);
    rho = reshape(rho,[1 2*ones(1,2*Q)]);
    
    trapped_ion_errormodel1_plus_idle(0)
    ECex = EC;
    trapped_ion_errormodel1_plus_idle(f)
    ECn = EC;

    for k = 1:size(algc,2)
        aa = algc{k};
        for n = 1:size(aa,1)
            
                
            
            if aa{n,3} ~= 'I'
                rho = apply_gate(rho,gate(aa{n,1},aa{n,3}),aa{n,2});
                
            end
            if any(ismember(aa{n,2},qnbs_n))
                rho = apply_error_channel(rho, ECn(['E' aa{n,3}]) , aa{n,2},Q);
            else
                rho = apply_error_channel(rho, ECex(['E' aa{n,3}]) , aa{n,2},Q);
            end
        end

    end

    rho = reshape(rho,2*ones(1,2*Q));
    %rho = reshape(rho,[1 2*ones(1,2*Q)]);

end
