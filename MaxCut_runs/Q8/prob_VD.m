% function to return the ideal probabilities neccesary to use VD
% these are done assuming perfect swaps, meaning one can directly
% computer rho^n
function [prob0_array,prob0p_array] = prob_VD(circ,max_copies,Q,f)
% circ is the circuit, n is the number of copies, Q is the number
% of qubits.
prob0_array = zeros(max_copies-1,1);
prob0p_array = zeros(max_copies-1,1);
for n = 2:max_copies
    rhon = frhon(circ,n,Q,f); % conpoute \rho^n (Q qubits required) 
    a = ftr_rhon_obs(rhon,Q);%compute Tr(\rho^n Z_1)
    b = ftr_rhon(rhon,Q);%compute Tr(\rho^n)
    % compute probabilit of getting 0  as a result of the ancilla measurement
    prob0_array(n-1) = (b+1)/2; %
    prob0p_array(n-1) = (a+1)/2;
end