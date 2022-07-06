function [Lalg,Lt] = to_layers(alg,Q)

    Lalg = {};
    Lt = zeros(size(alg,1),2);
    posto = 1:size(alg,1);

    ctrL = 0;
    while size(alg,1) > 0

        L = {}; %layer
        D = []; %description of a layer (qubits already taken)
        c = 0;
        ctr = 0;
        while c == 0
            c = 1;
            for q = 1:Q
                %find the first gate that touches q
                j = 1;
                while size(alg,1)>= j & ~any(alg{j,2} == q)
                    j = j+1;
                end
                %can I add that gate?
                if size(alg,1)>= j
                    if  size(alg{j,2},2) == 2
                        %two qubit gate
                        %are both qubits good to go?
                        c1 = 0;
                        for k = 1:j-1
                            if (any(alg{k,2} == alg{j,2}) == 1) | ...
                                (any(alg{k,2} == alg{j,2}(end:-1:1)) == 1)
                                c1 = 1;
                            end
                        end
                        if c1 == 0 & length(intersect(D,alg{j,2})) == 0
                            D = [D alg{j,2}];
                            c = 0;
                            L{end+1,1} = alg{j,1};
                            L{end,2} = alg{j,2};
                            L{end,3} = alg{j,3};
                            Lt(posto(j),:) = [ctrL+1,ctr+1];
                            %remove from alg
                            [alg,posto] = delete_from_alg(alg,j,posto);
                            ctr = ctr + 1;
                        end

                    else
                        %one qubit gate
                        if length(intersect(D,alg{j,2})) == 0
                            %add
                            D = [D alg{j,2}];
                            c = 0;
                            L{end+1,1} = alg{j,1};
                            L{end,2} = alg{j,2};
                            L{end,3} = alg{j,3};
                            Lt(posto(j),:) = [ctrL+1,ctr+1];
                            %remove from alg
                            [alg,posto] = delete_from_alg(alg,j,posto);
                            ctr = ctr + 1;
                        end
                    end
                end

            end
        end

        Lalg{end+1} = L;
        ctrL = ctrL+1;

    end


end
