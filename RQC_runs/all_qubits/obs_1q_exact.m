function [obs,rho] = obs_1q_exact(circ,Q)
    sigz = diag([1 -1]);
    circ = compile(circ,Q);
    f = 0;
    rho = init_rho(1:Q,[], f);
    rho = evol_alg(circ,rho,1:Q,[], f);
    rho = reshape(rho,[2 2^(Q-1) 2 2^(Q-1)]);
    obs = real(ncon({rho,sigz},{[2 1 3 1],[3 2]}));
    rho = reshape(rho,[2^Q 2^Q]);
end