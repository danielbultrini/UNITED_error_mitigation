path(pathdef)
clear all
seed=1;
nlsp=3;
f=1;
max_copies=6;
Q=6;
N=10;
Nts=100;
[circ,pos,ang] = maxcut_circ(seed,Q);
circ_inc_noise = increase_noise_level(circ, nlsp, seed);
graphdef = define_graph_rand(seed,Q);
qbsm  =   sort([graphdef(1,1), graphdef(1,2)],'ascend');
circnC = gen_ts(circ,pos,ang,N,Nts,seed);
[obsex_coi,rhoexact] = obs_exact(circ,Q,qbsm);
 [obsn_coi, prob0_array_coi,...
            prob0p_array_coi,rhonoisy,circc_coi] = obs(circ_inc_noise, Q, max_copies,f,qbsm);
circnC = gen_ts(circ,pos,ang,N,Nts,seed);
for c = 1:length(circnC)
        tcirc = circnC{c};
        % increase the noise:
        tcirc_inc_noise = increase_noise_level(tcirc, nlsp, seed);
        % calculate the useful information:
            [obsn_train, prob0_train,...
        prob0p_train,~,circc]= obs(tcirc_inc_noise, Q, max_copies,f, qbsm);
        obsex = obs_exact(tcirc,Q, qbsm);
        % store the data:
        obsn_array_train( c) = obsn_train;
        obsex_array_train( c) = obsex;
        prob0_array_train(c,:) = prob0_train;
        prob0p_array_train(c,:) = prob0p_train;
        circc_train{c} = circc;
        compt(c) = compare(circc_train{c},circc_coi);
end
min(compt)

figure(); hold on;
plot(obsex_array_train,obsn_array_train,'*')
plot(obsex_array_train,(2*prob0p_array_train(:,1)-1)./(2*prob0_array_train(:,1)-1),'+')
plot(obsex_array_train,(2*prob0p_array_train(:,5)-1)./(2*prob0_array_train(:,5)-1),'d')
plot(obsex_coi,obsn_coi,'r*')
plot(obsex_coi,(2*prob0p_array_coi(1)-1)./(2*prob0_array_coi(1)-1),'r+')
plot(obsex_coi,(2*prob0p_array_coi(5)-1)./(2*prob0_array_coi(5)-1),'rd')


function comp = compare(circ1,circ2)
  comp=1;
  if length(circ1) == length(circ2)
    for i = 1:length(circ1)
      if size(circ1{i},1) == size(circ2{i},1)
         for ii = 1:size(circ1{i},1)
           if ~isequal(circ1{i}{ii,2},circ2{i}{ii,2})
             comp =0;
           end
           if ~isequal(circ1{i}{ii,3},circ2{i}{ii,3})
             comp =0;
           end
         end
      else
        comp=0;     
      end
    end
  else
   comp=0;  
  end
end
