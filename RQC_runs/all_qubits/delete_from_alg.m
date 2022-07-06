function [algt,posto] = delete_from_alg(alg,p,posto)
    algt = {};
    for i = 1:p-1
        for j = 1:3
            algt{i,j} = alg{i,j};
        end
    end

    for i = p+1:size(alg,1)
        for j = 1:3
            algt{i-1,j} = alg{i,j};
        end

    end

   posto = posto([1:p-1 p+1:size(alg,1)]);     

end
