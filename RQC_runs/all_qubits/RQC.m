function [alg,pos,angles] = RQC(p,Q,seed)

rng(seed)
qnbs = 1:Q;
% creates RQC with p layers for Q qubits
% the layout from Fig. 3 of REQUEST paper  
% parameters values determined by seed

% pos a structure carrying information about  types of gates in alg 
% angles a structure carrying information about angles of gates in alg 

alg = {};

Qrange1 = 1:2:Q-1;
Qrange1s1 = min(Qrange1):max(Qrange1)+1;
l1g2 = length(Qrange1);
l1g1 = length(Qrange1s1);

Qrange2 = 2:2:Q-1;
Qrange2s1 = min(Qrange2):max(Qrange2)+1;
l2g2 = length(Qrange2);
l2g1 = length(Qrange2s1);

pos = struct;
angles = struct;

pos.x = [];
pos.y = [];
pos.z = [];
pos.xx  = [];

angles.x = [];
angles.y = [];
angles.z = [];
angles.xx  = [];

ctr = 0;

for ii = 1:p
    ang = rand(l1g1,3)*2*pi;
    qnbsl = qnbs(Qrange1s1);
    for i = 1:l1g1
        [algt,pos,angles,ctr] = ...
            u3(qnbsl(i),ang(i,1),ang(i,2),ang(i,3),pos,angles,ctr);
        alg = vertcat(alg,algt);
    end
    
    ang = rand(l1g2,1)*pi;
    qnbsl = qnbs(Qrange1);
    for i = 1:l1g2
        alg = vertcat(alg,{ang(i),[qnbsl(i) qnbsl(i)+1],'2'});
        ctr = ctr+1;
        pos.xx(end+1) = ctr;
        angles.xx(end+1) = ang(i);
    end
    
    ang = rand(l2g1,3)*2*pi;
    qnbsl = qnbs(Qrange2s1);
    for i = 1:l2g1
        [algt,pos,angles,ctr] = ...
            u3(qnbsl(i),ang(i,1),ang(i,2),ang(i,3),pos,angles,ctr);
        alg = vertcat(alg,algt);
    end
    
    ang = rand(l2g2,1)*pi;
    qnbsl = qnbs(Qrange2);
    for i = 1:l2g2
        alg = vertcat(alg,{ang(i),[qnbsl(i) qnbsl(i)+1],'2'});
        ctr = ctr+1;
        pos.xx(end+1) = ctr;
        angles.xx(end+1) = ang(i);
    end   
end

end

function [algc, pos, ang, ctr] = u3(qnb,theta,phi,lam, pos, ang, ctr)

% replace 1 qubit gates by the native gates
% U(\theta,\phi,\lambda) = [ [cos(\theta/2) -sin(\theta/2)e^{i\lambda}];
%           [e^{i\phi} sin(\theta/2)  e^{i(\phi+\lambda)} sin(\theta/2)]]; 
% U(\theta,\phi,\lambda) = R_Z(\phi) R_Y(\theta) R_Z(\lambda)
% R_Z(\theta) = e^{1i \theta/2 \sigma_Z}
% decompose general single unitary 

algc = {};
algc =  vertcat(algc,{lam/2,qnb,'Z'});
algc =  vertcat(algc,{theta/2,qnb,'Y'});
algc =  vertcat(algc,{phi/2,qnb,'Z'});

pos.z(end+1) = ctr+1;
pos.y(end+1) = ctr+2;
pos.z(end+1) = ctr+3;

ang.z(end+1) = lam/2;
ang.y(end+1) = theta/2;
ang.z(end+1) = phi/2;

ctr = ctr+3;

end