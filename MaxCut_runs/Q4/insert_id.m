function Lalg2 = insert_id(Lalg,Q)


    Lalg2 = Lalg;  




    for k = 1:size(Lalg,2)
        aa = Lalg{k};
       

        %are there "empty" qubits in a layer?
        %if yes, apply idle error channel
        D = [];
        for n = 1:size(aa,1)
            D = [D aa{n,2}];
        end
        D = setdiff([1:Q],D);
        
        for n = D
           Lalg2{k} = vertcat(Lalg2{k},{'-',n,'I'});
        end

    end






end