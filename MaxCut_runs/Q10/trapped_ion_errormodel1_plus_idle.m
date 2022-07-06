function trapped_ion_errormodel1_plus_idle(f)
    global EC ECrates

    % paulis
    I = eye(2);
    X = [0 1; 1 0];
    Y = [0 -1.i; 1.i 0];
    Z = [1 0; 0 -1];


   
    % realistic error rates
    p_alpha = 1e-4*f ;
    p_dep = 8e-4*f   ;
    p_xx = 1e-3*f    ;
    p_h = 1.25e-3*f  ;
    p_d = 1.5e-4*f   ;
    p_d2 = 7.5e-4*f  ;
    p_idle = 8e-4*f  ;


    [EX,errX] = nX(0, p_d, p_dep, p_alpha);EX = reshape(EX,2,2,2,2);
    [EY,errY] = nY(0, p_d, p_dep, p_alpha);EY = reshape(EY,2,2,2,2);
    [EZ,errZ] = nZ(0, p_d, p_dep, p_alpha);EZ = reshape(EZ,2,2,2,2);
    [EI,errI] = nIdle(p_idle);EI = reshape(EI,2,2,2,2);
    [E2,err2] = nXX(0, p_d2, p_dep, p_h, p_xx);
    E2 = reshape(E2,2,2,2,2,2,2,2,2);
    E2 = permute(E2,[1 3 2 4 5 7 6 8]);

    EC = containers.Map;
    EC('EX') = EX;
    EC('EY') = EY;
    EC('EZ') = EZ;
    EC('EI') = EI;
    EC('E2') = E2;
    
    ECrates = containers.Map;
    ECrates('EX') = errX;
    ECrates('EY') = errY;
    ECrates('EZ') = errZ;
    ECrates('EI') = errI;
    ECrates('E2') = err2;    






    %% define native state preparation
    % ground state preparation for a single qubit
    gs = [1 0 0 0]';    % ideal ground state (vectorized density matrix)
    ngs = singleq_depol_channel(p_dep)*gs;      % depolarized ground state
    EC('ngs') = ngs;
    ECrates('ngs') =  p_dep;

    %% define native POVM elements
    % vectorized POVM effect for measurement of ground state
    E0 = [1 0; 0 0];    % ideal ground state projector
    noisy_E0 = singleq_depol_channel(p_dep) * reshape(E0, [4,1]); % apply depolarizing channel
    nE0 = noisy_E0';        % Take its adjoint.
    EC('nE0') = nE0;
    ECrates('nE0') =  p_dep;
    

    % Why does this look like this?
    % The error model says that the noisy measurement is the ideal measurement
    % preceded by a depolarizing channel. So the outcome probabilities for a
    % state \rho look like:
    %   tr( E D(\rho) E'), where D(.) is the depolarizing channel, and E are
    %     the ideal POVM effects.
    % This becomes:
    %   tr( E'E D(\rho) ) = tr( D(E'E) \rho) since the depolarizing map is it's
    %   own adjoint.
    % Writing this probability in terms of vectorized quantities gives us:
    %   tr( D(E'E) \rho) = ( D(E'E) | \rho ),
    %     where |o) is the vectorized version of the matrix o and (o| is it's
    %     adjoint. So we need ( D(E'E)| to extract probabilites, which is what
    %     the above calculates

    % vectorized POVM effect for measurement of excited state
    E1 = [0 0; 0 1];    % ideal excited state projector
    noisy_E1 = singleq_depol_channel(p_dep) * reshape(E1, [4,1]); % apply depolarizing channel
    nE1 = noisy_E1';        % Take its adjoint.
    EC('nE1') = nE1;
    ECrates('nE1') =  p_dep;


    %% define native one qubit gates -- these are parameterized by the ideal
    %% rotation angle. The functions below define the map for any ideal rotation, and then
    %% we instantiate commonly used gates, i.e.,
    %% nX(pi/2,...), nX(-pi/2,...), nX(pi,...)
    %% nY(pi/2,...), nY(-pi/2,...), nY(pi,...)
    %% nZ(pi/2,...), nZ(pi,...)

    nX_pi2 = nX(pi/2, p_d, p_dep, p_alpha);
    nX_mpi2 = nX(-pi/2, p_d, p_dep, p_alpha);
    nX_pi = nX(pi, p_d, p_dep, p_alpha);
    nY_pi2 = nY(pi/2, p_d, p_dep, p_alpha);
    nY_mpi2 = nY(-pi/2, p_d, p_dep, p_alpha);
    nY_pi = nY(pi, p_d, p_dep, p_alpha);
    nZ_pi2 = nZ(pi/2, p_d, p_dep, p_alpha);
    nZ_pi = nZ(pi, p_d, p_dep, p_alpha);


    % native idle
    function [G,p_dep] = nIdle(p_dep)
        I = eye(2);
        G = singleq_depol_channel(p_dep) * ...                         % depolarizing channel
        kron(I, I) ;                                             % the ideal gate
    end

    % native X rotation
    function [G,p] = nX(theta, p_d, p_dep, p_alpha)
        % note: this theta is the QI convention for rotation angle -- i.e., theta/2
        % is the angle that is exponentiated.
        X = [0 1; 1 0];
        G = singleq_Z_channel(p_d) * ...                               % dephasing channel
        singleq_depol_channel(p_dep) * ...                         % depolarizing channel
        singleq_X_channel(p_alpha) * ...                           % imprecise rotation
        kron(transpose(expm(1.i*(theta/2)*X)), expm(-1.i*(theta/2)*X)) ; % the ideal gate
        p= p_d+p_dep+p_alpha;
    end

    % native Y rotation
    function [G,p] = nY(theta, p_d, p_dep, p_alpha)
        % note: this theta is the QI convention for rotation angle -- i.e., theta/2
        % is the angle that is exponentiated.
        Y = [0 -1.i; 1.i 0];
        G = singleq_Z_channel(p_d) * ...                               % dephasing channel
        singleq_depol_channel(p_dep) * ...                         % depolarizing channel
        singleq_Y_channel(p_alpha) * ...                           % imprecise rotation
        kron(transpose(expm(1.i*(theta/2)*Y)), expm(-1.i*(theta/2)*Y)) ; % the ideal gate
        p = p_d+p_dep+p_alpha;
        
    end

    % native Z rotation
    function [G,p] = nZ(theta, p_d, p_dep, p_alpha)
        % note: this theta is the QI convention for rotation angle -- i.e., theta/2
        % is the angle that is exponentiated.
        Z = [1 0; 0 -1];
        G = singleq_Z_channel(p_d) * ...                               % dephasing channel
        singleq_depol_channel(p_dep) * ...                         % depolarizing channel
        singleq_Z_channel(p_alpha) * ...                           % imprecise rotation
        kron(transpose(expm(1.i*(theta/2)*Z)), expm(-1.i*(theta/2)*Z)) ; % the ideal gate
         p = p_d+p_dep+p_alpha;
    end

    %% define native two qubit gate
    function [G,p] = nXX(theta, p_d2, p_dep, p_h, p_xx)
        X = [0 1;1 0];
        XX = kron(X,X);
        G = kron(singleq_Z_channel(p_d2), singleq_Z_channel(p_d2)) * ...             % dephasing channel
        kron(singleq_depol_channel(p_dep), singleq_depol_channel(p_dep)) * ... % depolarizing channel
        twoq_XX_channel(p_h) * ...                                             % heating error channel
        twoq_XX_channel(p_xx) * ...                                            % imprecise rotation
        kron(transpose(expm(1.i*(theta/2)*XX)) , expm(-1.i*(theta/2)*XX));     % the ideal gate
        p = 2*p_d2+2*p_dep+p_h+p_xx;
    end


    %% define error channels
    function E = singleq_X_channel(p)
        I = eye(2); X = [0 1; 1 0];

        E = (1-p)*kron(I,I) + (p)*kron(X,X);
    end

    function E = singleq_Y_channel(p)
        I = eye(2); Y = [0 -1.i; 1.i 0];

        E = (1-p)*kron(I,I) + (p)*kron(transpose(Y),Y);
    end

    function E = singleq_Z_channel(p)
        I = eye(2); Z = [1 0; 0 -1];

        E = (1-p)*kron(I,I) + (p)*kron(Z,Z);
    end

    function E = singleq_depol_channel(p)
        I = eye(2); X = [0 1; 1 0]; Y = [0 -1.i; 1.i 0]; Z = [1 0; 0 -1];

        E = (1-p)*kron(I,I) + (p/3)*kron(X,X) + (p/3)*kron(transpose(Y),Y) + (p/3)*kron(Z,Z);
    end

    function E = twoq_XX_channel(p)
        I = eye(2); X = [0 1; 1 0];
        XX = kron(X,X);
        II = kron(I,I);

        E = (1-p)*kron(II,II) + (p)*kron(transpose(XX),XX);
    end
end
