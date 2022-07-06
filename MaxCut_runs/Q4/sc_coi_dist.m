clear all
Q=8;
for seed = 1:30
try 
[circ,pos,ang,qbsm] = maxcut_circ(seed,Q);
[obs(seed)] = obs_exact(circ,Q,qbsm);
catch 
disp(seed)
end
end
sort(abs(obs),'ascend')

