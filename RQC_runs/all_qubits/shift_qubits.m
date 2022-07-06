function algn = shift_qubits(alg,Qshift)
    % shfts a circuit by Q shift qubits 
    % e. g. for Qshift = 2 a circuit acting at qubits 1..4 is replaced by
    % a circuit acting at qubits 3..6
    
    algn = alg;
    for i = 1:size(alg,1)
        algn{i,2} = alg{i,2} + Qshift;
    end   

end