function a = ftr_rhon(rho,Q)
    a = real(trace(reshape(rho,[2^Q 2^Q])));
end