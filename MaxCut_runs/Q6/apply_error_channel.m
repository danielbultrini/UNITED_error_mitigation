function y = apply_error_channel(x,g,p,Q)


    sv = size(x);
        

    if length(p) == 1
        % 1 qubit        
        p1  = [0 1:(p-1) (p+1):(p+Q-1) (p+Q+1):(2*Q) p p+Q]+1;
        x = reshape(permute(x,p1),[2^(2*Q-2) 4]);
        y = reshape(x*reshape(g,[4 4]),sv);
        p1 = [0 1:(p-1) 2*Q-1 p:(p+Q-2) 2*Q (p+Q-1):(2*Q-2)]+1;
        y = permute(y,p1); 
    else
        if p(2) < p(1)
            p = p([2 1]);
            g = permute(g,[2 1 4 3 6 5 8 7]);
        end
        % 2 qubit
        p1  = [0 1:(p(1)-1) p(1)+1:(p(2)-1) (p(2)+1):(p(1)+Q-1) ...
            (p(1)+Q+1):(p(2)+Q-1) (p(2)+Q+1):(2*Q) p p+Q]+1;
        x = reshape(permute(x,p1),[2^(2*Q-4) 16]);
        y = reshape(x*reshape(g,[16 16]),sv);
         p1  = [0 1:(p(1)-1) 2*Q-3 p(1):(p(2)-2) 2*Q-2 (p(2)-1):(p(1)+Q-3) ...
            2*Q-1 (p(1)+Q-2):(p(2)+Q-4) 2*Q (p(2)+Q-3):(2*Q-4)]+1;
        y = permute(y,p1); 
    end



end
