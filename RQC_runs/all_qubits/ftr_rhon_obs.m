function a = ftr_rhon_obs(rho,Q)
    sigz = diag([1 -1]);
    rho = reshape(rho,[2 2^(Q-1) 2 2^(Q-1)]);
    a = real(ncon({rho,sigz},{[2 1 3 1],[3 2]}));
end