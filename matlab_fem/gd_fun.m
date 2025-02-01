function [valore] = gd_fun(x,y,t, marker)
        
        gd_fun1 = @(x,y,t) sin(t)+y^2+1;
        gd_fun3 = @(x,y,t) sin(t)+x^2+1;
        gd_fun5 = @(x,y,t) sin(t)+y^2;
        
        if marker == 1
            valore = gd_fun1(x,y,t);
        elseif marker == 3 
            valore = gd_fun3(x,y,t);
        else
            valore = gd_fun5(x,y,t);
        end
end