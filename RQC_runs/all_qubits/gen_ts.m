function [algproj] = gen_ts(alg,pos,ang,N,Nts,seed)

% generates training set for a circuit alg with gates positions pos and
% angles ang
% algproj{i} is a cell array with Nts near-Clifford circuits
% each of them has N non Cliffor gates 
% seed controls the projection



rng(seed)
sigma = 0.5;

for i = 1:Nts
    angproj{i} = proj_angles_all(ang,sigma,N);
    algproj{i} = angles_to_alg(alg, pos, angproj{i});
end

end


function algproj = angles_to_alg(alg, pos, ang)

algproj = alg;

for i = 1:length(pos.x)
    algproj{pos.x(i),1} = ang.x(i);
end

for i = 1:length(pos.y)
    algproj{pos.y(i),1} = ang.y(i);
end

for i = 1:length(pos.z)
    algproj{pos.z(i),1} = ang.z(i);
end

for i = 1:length(pos.xx)
    algproj{pos.xx(i),1} = ang.xx(i);
end


end


function angproj = proj_angles_all(ang,sigma,nbnonCliff)

    angproj = struct;
    angproj.x  =  project_angles('X',ang.x,sigma,0);
    angproj.y  =  project_angles('Y',ang.y,sigma,0);
    angproj.z  =  project_angles('Z',ang.z,sigma,nbnonCliff);
    angproj.xx  =  project_angles('2',ang.xx,sigma,0);

end


function  [alphaproj] = project_angles(fl,alpha,sigma,nbnonCliff)

% Projects a vector of angles of the fl = Z,X,Y,XX rotations alpha, e. g.
% R_{fl}(\alpha(i))  = expm(-1i*alpha(i)*fl),
%(a non standard R_{fl} rotaton definition used here !!!).
% alpha is consistent with the angles apperaing in Lukasz's circuits
% parametrization.

% Return angles alphaproj of the projected gates.
% nbnonCliff is a number of angles (gates) left intact by the projection.
% A likelihood of a gate being projected to a  Clifford rotation R_{fl}(k*pi/4), k = 0,1,2,3  
% is proprtional to exp(- d/\sigma^2), 
%d = ||e^{i\alpha}R_{fl}(alpha) - e^{ik\pi/4}R_{fl}(k*pi/4)||/d, 
% with ||.|| being  Frobenius norm. 
% here || R_Z(aplha)/d|| = 1

    for ii = 1:4
        gateref{ii} = gate_gauge((ii-1)*pi/4,fl);
    end
    
    pflipgt = [];
    pflipgc = {};
    ll = length(alpha);
    alphaproj = alpha;
    for i = 1:ll
                gate1 = gate_gauge(alpha(i),fl);          
                for  ii = 1:4          
                    dist(ii) = norm(gate1-gateref{ii});
                end
                prob = exp(-(dist/sigma).^2);
                pflipgc{i} = prob;
                pflipgt(i) = sum(pflipgc{i});
    end
    flipgt = [];
    for i = 1:ll-nbnonCliff
        cum = cumsum(pflipgt/sum(pflipgt));
        x = rand();
        outcome = x>  cum;
        outcome = sum(outcome)+1;
        flipgt(i) = outcome;
        pflipgt(outcome) = 0;
    end
    for i = 1:ll-nbnonCliff
        x = rand();
        ii = flipgt(i);
        cum = cumsum(pflipgc{ii}/sum(pflipgc{ii}));
        outcome = x>  cum;
        outcome = sum(outcome);
        ang = (outcome)*pi/4;
        alphaproj(ii) = ang;
    end
    
end


function m = gate_gauge(a,fl)

    if fl == 'X'
        m = cos(a)*eye(2) - sin(a)*[0 1;1 0]*1.i;
        m = exp(1i*a)*m;
        m = reshape(m, [4 1])/sqrt(2);
    elseif fl == 'Y'
        m = cos(a)*eye(2) - sin(a)*[0 -1.i;1.i 0]*1.i;
        m = exp(1i*a)*m;
        m = reshape(m, [4 1])/sqrt(2);
    elseif fl == 'Z'
        m = cos(a)*eye(2) - sin(a)*[1 0;0 -1]*1.i;
        m = exp(1i*a)*m;
        m = reshape(m, [4 1])/sqrt(2);
    elseif fl == '2'
        m = cos(a)*eye(4) - sin(a)*kron([0 1;1 0],[0 1;1 0])*1.i;
        m = exp(1i*a)*m;
        m = reshape(m, [16 1])/2;
    end

end


