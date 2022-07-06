function [obs,rho] = obs_exact(circ,Q,ii)
    sigx = [[0 1]; [1 0]];
    circ = compile(circ,Q);
    f = 0;
    rho = init_rho(1:Q,[], f);
    rho = evol_alg(circ,rho,1:Q,[], f);
    obs = ncon({rho,sigx,sigx},...
   {[1:ii(1)-1 Q+1 ii(1)+1:ii(2)-1 Q+3 ii(2)+1:Q 1:ii(1)-1 Q+2 ii(1)+1:ii(2)-1 Q+4 ii(2)+1:Q],[Q+1 Q+2],[Q+3 Q+4]});
    obs = real(obs);
end
