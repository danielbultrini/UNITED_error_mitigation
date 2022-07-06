function y = apply_right(x,g,p)
    %applies gate g at positions p "to the right" (normal direction)
    %x represents input states: first index enumerates states. After that,
    %the indices are numbered from top to bottom.
    %index convention for g: first left, from top

    N = length(size(x)); %number of indices of x
    q = length(p); %number of qubits that g acts on
    q2 = 2^q;


    %build the first permutation
    t = 1:N;
    t = t(~ismember(1:N,p+1));
    per = [t p+1];
    y = permute(x, per);
    y = reshape(y,numel(y)/q2,q2);
    y = y * reshape(g,q2,q2);
    y = reshape(y,[size(x,1) 2*ones(1,N-1)]);
    %the second permutation is the inverse of the first one
    per2(per) = 1:length(per);
    y = permute(y,per2);


end
