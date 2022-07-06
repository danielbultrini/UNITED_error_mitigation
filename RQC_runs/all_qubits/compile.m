function algc = compile(alg,Q)
  % returns a compiled circuit with layered strucure and identities inserted 
  
  [algc,~] =   to_layers(alg,Q);
  algc = insert_id(algc,Q);
  
  
end