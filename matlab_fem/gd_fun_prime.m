function [valore] = gd_fun_prime(x,y, t,  marker)

        
        gd_fun1 = @(x,y, t) cos(t);
        gd_fun3 = @(x,y, t) cos(t);
        gd_fun5 = @(x,y, t) cos(t);
        
        if marker == 1
            valore = gd_fun1(x,y,t);
        elseif marker == 3 
            valore = gd_fun3(x,y,t);
        else
            valore = gd_fun5(x,y,t);
        end
end