function m = gate(a,fl)

    if fl == 'X'
        m = cos(a)*eye(2) - sin(a)*[0 1;1 0]*1.i;
    elseif fl == 'Y'
        m = cos(a)*eye(2) - sin(a)*[0 -1.i;1.i 0]*1.i;
    elseif fl == 'Z'
        m = cos(a)*eye(2) - sin(a)*[1 0;0 -1]*1.i;
    elseif fl == 'I'
        m = eye(2);
    elseif fl == '2'
        m = cos(a)*eye(4) - sin(a)*kron([0 1;1 0],[0 1;1 0])*1.i;
        m = reshape(m,2,2,2,2);
    elseif fl == 'ga'
        m = a;
    end

end
