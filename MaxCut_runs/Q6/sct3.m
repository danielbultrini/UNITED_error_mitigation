path(pathdef)
% layer number 
p = 4;
% number of qubits 
Q=6;
% seed 
seed = 1;
% random quantum cicuit
[circ,~,~] = RQC(p,Q,seed);
% number of copies
n = 2;
% CTRL-SWAP noise level
nl = 1;

% observable Z_1
% nl = 1 corresponds to diagrams from Fig. 2
% nl = 3, 5, 7, 9 ... perform identity insertion 

%compute Tr(\rho^n Z_1)/Tr(\rho^n) using a circuit with perfect swaps  (2*Q+1 qubits required) 
% [prob0,prob0p] = VD_ideal(circ,n,Q,nl);
% obsm = (2*prob0p-1)/(2*prob0-1); % 
f = 1/length(circ);
rhon = frhon(circ,n,Q,f); % conpoute \rho^n (Q qubits required) 
a = ftr_rhon_obs(rhon,Q);%compute Tr(\rho^n Z_1)
b = ftr_rhon(rhon,Q);%compute Tr(\rho^n)
% compute probabilit of getting 0  as a result of the ancilla measurement
prob0v2 = (b+1)/2; %
prob0pv2 = (a+1)/2;
obsm = (2*prob0pv2-1)/(2*prob0v2-1);
% check consistency  of both ways of computing the probabilities


% % compute the exact and noisy expectation values for single copy 
[obsex,rhoex] = obs_1q_exact(circ,Q);
[obsn,rhon] = obs_1q_noisy(circ,Q,f);

%[obsex, obsm, obsn, obsm-obsm2] 

 

    

 
    
