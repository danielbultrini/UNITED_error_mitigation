% function to calculate the noisy observable
function [obsn, prob0_array, prob0p_array,rhon] = obs(circ, Q, max_copies, f)

    rhon = frhon_eff(circ,max_copies,Q,f); % conpoute \rho^n (Q qubits required) 
    obsn = obs_1q_noisy_eff(rhon{1},Q);
    
    [prob0_array, prob0p_array] = prob_VD_eff(rhon,max_copies,Q);
end

function [prob0_array,prob0p_array] = prob_VD_eff(rhon,max_copies,Q)
% circ is the circuit, n is the number of copies, Q is the number
% of qubits.
prob0_array = zeros(max_copies-1,1);
prob0p_array = zeros(max_copies-1,1);
for n = 2:max_copies
   
    a = ftr_rhon_obs(rhon{n},Q);%compute Tr(\rho^n Z_1)
    b = ftr_rhon(rhon{n},Q);%compute Tr(\rho^n)
    % compute probabilit of getting 0  as a result of the ancilla measurement
    prob0_array(n-1) = (b+1)/2; %
    prob0p_array(n-1) = (a+1)/2;
end

end

function [obs, rho] = obs_1q_noisy_eff(rho,Q)
    sigz = diag([1 -1]);
    rho = reshape(rho,[2 2^(Q-1) 2 2^(Q-1)]);
    obs = real(ncon({rho,sigz},{[2 1 3 1],[3 2]}));
    rho = reshape(rho,[2^Q 2^Q]);
end

function [rhon]  = frhon_eff(circ,n_max,Q,f)

% compute rho^n for noisy rho defined by circuit rho 
% Q is number of the qubits 

circc = compile(circ,Q);
rho = init_rho([],1:Q, f);
rho = evol_alg(circc,rho,[],1:Q, f);
rho = reshape(rho,[2^Q 2^Q]); 
rhon{1} = rho;
for i = 2:n_max
    rhon{i} = rhon{i-1}*rho;
end
for i = 1:n_max
rhon{i} = reshape(rhon{i},2*ones(1,2*Q));
end
end
