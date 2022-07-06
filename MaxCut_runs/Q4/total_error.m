function errtot = total_error(circ,Q)

global ECrates 

trapped_ion_errormodel1_plus_idle(1)
circc = compile(circ,Q);
errtot = Q*ECrates('ngs');


for i = 1:length(circc)
    for j = 1:size(circc{i},1)
         errtot = errtot + ECrates(['E' circc{i}{j,3}]);
    end
end


end