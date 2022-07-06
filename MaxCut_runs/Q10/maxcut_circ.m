function [alg,pos,angles,qbsm] = maxcut_circ(seedg,Q) 
% QAOA for a maxcut random graph Hamiltonian as described in 
% NIBP paper
% the cost function is EV of H defined in the paper given by mxH
% the function is suppposed to be efficient for small number of qubits n 
% for which scon overhead is dominant 
%  |psi> = ...  e^(-i param(1) H_X) e^(-i param(2) H_ZZ)  |+>
% ord has the samme order as param
% line connectivity (1,2),(2,3),...,(Q-1,Q)
graphdef = define_graph_rand(seedg,Q);
in=load(sprintf('angles_Q%d_seedg%d',Q,seedg));
param=in.angles;
qbsm=in.qbsm;

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

p =  length(param)/2;
alg = {};
for i = 1:p
    [algt,pos,angles,ctr] =  lay_qaoa(param(2*(i-1)+1),param(2*i),Q,graphdef,pos,angles,ctr);
    alg  = vertcat(alg, algt);  
end
      
end

function [alg,pos,angles,ctr]  = lay_qaoa(gamma,beta,Q,graphdef,pos,angles,ctr)

alg = {};

for i = 1:size(graphdef,1)
    alg  = vertcat(alg,{gamma,[graphdef(i,1),graphdef(i,2)],'2'});
    ctr = ctr +1;
    pos.xx(end+1) = ctr;
    angles.xx(end+1) = gamma;
end

for i = 1:Q
    alg  = vertcat(alg,{beta,i,'Z'});
    ctr = ctr+1;
    pos.z(end+1) = ctr;
    angles.z(end+1) = beta;
end

end
