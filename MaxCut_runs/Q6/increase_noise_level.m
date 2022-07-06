function new_circ = increase_noise_level(circ, nl, seed)
% function to increase noise in circuit by indentity, Clifford
%insertions for each gate. 
    rng(seed)
    rand_ints = randi([0,3], length(circ));
    if nl ==0
        new_circ = circ;
    elseif nl ==1
        new_circ = circ;
    else
        new_circ = {};
        for i=1:length(circ)
            ang = circ{i,1};
            qubits = circ{i,2};
            gate = circ{i,3};
            % check if 2 qubit gate for rand clifford.
            if gate == '2'
                rand_cliff = pi/4*rand_ints(i);
            else
                rand_cliff = pi/2*rand_ints(i);
            end
            for n=1:nl
                if n ~= nl
                    new_lay = {(-1)^(n-1)*rand_cliff, qubits, gate};
                    new_circ = [new_circ; new_lay];
                else 
                    if rem(nl,2)==0 
                        new_lay = {ang-rand_cliff, qubits, gate};
                        new_circ = [new_circ; new_lay]; 
                    else
                        new_lay = {ang, qubits, gate};
                        new_circ = [new_circ; new_lay]; 
                    end
                end    
            end
        end
    end
end

