function obs = ftr_rhon_obs(rho,Q,ii)
    sigx = [[0 1]; [1 0]];
    obs = ncon({rho,sigx,sigx},...
   {[1:ii(1)-1 Q+1 ii(1)+1:ii(2)-1 Q+3 ii(2)+1:Q 1:ii(1)-1 Q+2 ii(1)+1:ii(2)-1 Q+4 ii(2)+1:Q],[Q+1 Q+2],[Q+3 Q+4]});
    obs = real(obs);
end
