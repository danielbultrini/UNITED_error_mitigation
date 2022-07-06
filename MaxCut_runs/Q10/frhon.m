function [rhon]  = frhon(circ,n,Q,f)

% compute rho^n for noisy rho defined by circuit rho 
% Q is number of the qubits 

circc = compile(circ,Q);
rho = init_rho([],1:Q, f);
rho = evol_alg(circc,rho,[],1:Q, f);
rho = reshape(rho,[2^Q 2^Q]); 
rhon = rho;
for i = 2:n
    rhon = rhon*rho;
end
rhon = reshape(rhon,2*ones(1,2*Q));
end

