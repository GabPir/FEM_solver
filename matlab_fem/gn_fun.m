function [valore] = gn_fun(x,y,t, marker)

    gn_fun2 = @(x,y, t) -2*y;

    valore = gn_fun2(x,y,t);

end

