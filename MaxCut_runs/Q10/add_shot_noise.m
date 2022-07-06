function exp_w_shot_noise = add_shot_noise(exp, shots)
    p = (exp+1)/2;
    smax = 10^7; % size of the largest allocated array
    ss = 0;
    nl = 0;
    while ss < shots
        stmp = min(shots-ss,smax);
        rand_numbers = rand(stmp,1);
        mask = rand_numbers(rand_numbers < p);
        nl = nl + length(mask);
        ss = ss+stmp;
    end
    exp_w_shot_noise = 2*nl/shots - 1;
end

