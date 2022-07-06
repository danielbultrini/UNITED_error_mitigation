function y = apply_gate(x,g,p)


    N = size(size(x),2);

    y = apply_right(x,g,p + (N-1)/2 );
    y = apply_left(y,dagger(g),p);




end
