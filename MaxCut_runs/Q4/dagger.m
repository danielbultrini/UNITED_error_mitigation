function v = dagger(u)
    a = length(size(u));
    per = [a/2+1:a 1:a/2];
    v = permute(conj(u),per);
end
