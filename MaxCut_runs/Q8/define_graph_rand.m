function graphdef = define_graph_rand(seed,Q)

rng(seed)

ctre = 0;
for i = 1:Q
    for j = i+1:Q
        ctrl = round(rand());
        if ctrl 
            ctre = ctre +1;
            graphdef(ctre,:) = [i j]; 
        end
    end
end

end

