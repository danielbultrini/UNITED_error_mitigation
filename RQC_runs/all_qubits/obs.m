% function to calculate the noisy observable
function [obsn, prob0_array, prob0p_array,rho] = obs(circ, Q, max_copies, f)
    [obsn,rho] = obs_1q_noisy(circ,Q,f);  
    [prob0_array, prob0p_array] = prob_VD(circ,max_copies,Q, f);
end