function check_if_Cliff(alg)
    nCliff = 0;
    nnCliffrot = 0;
    nnCliff = 0;
%    angles_all = [];


   
    for i = 1:size(alg,1)
        if strcmp(alg{i,3},'X')  || strcmp(alg{i,3},'Y') ...
            || strcmp(alg{i,3},'Z') || strcmp(alg{i,3},'2')  
               if  ~(mod(4*alg{i,1},pi) < 1.e-10 || mod(4*alg{i,1},pi) > pi - 1.e-10)
                    nnCliffrot = nnCliffrot+1;
               else
                    nCliff = nCliff + 1;
               end

              
        else
            nnCliff =  nnCliff +1;
        end
    end
    str = ['number of Cliffords: ',num2str(nCliff),' number of non-Clifford rotations: ',num2str(nnCliffrot),...
        '  number of other gates: ',num2str(nnCliff)];
    disp(str);
end