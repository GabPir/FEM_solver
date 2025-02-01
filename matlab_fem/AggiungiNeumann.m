function [] = AggiungiNeumann()
    global geom gn pivotMidNe pivotMid ni t;
    
    % pesi quadratura nel segmento -1,1 trasformati nel caso 0,1
    % (accuratezza fino ad ordine 5)

    t_hat = [0.046910077030668, 0.230765344947158, 0.5, 0.769234655052842, 0.953089922969332];
    w = [0.118463442528095, 0.239314335249683, 0.284444444444444, 0.239314335249683, 0.118463442528095];



    % Valutazione delle funzioni nei possibili punti
    phib = 2*t_hat(1,:).^2 - 3.*t_hat(1,:) +1;  %phi beginning
    phim = -4.*t_hat(1,:).^2 + 4.*t_hat(1,:); %phi midpoint
    phie = 2.*t_hat(1,:).^2 - t_hat(1,:); %phi end
    




    Nb = length(geom.pivot.Ne(:,1));
    for b =1:Nb
        Vb = geom.elements.borders(geom.pivot.Ne(b,1),1); 
        Ve = geom.elements.borders(geom.pivot.Ne(b,1),2);
        E_norm = norm(geom.elements.coordinates(Vb,:)-geom.elements.coordinates(Ve,:),2);
        
        xb = geom.elements.coordinates(Vb,1);
        yb = geom.elements.coordinates(Vb,2);
        xe = geom.elements.coordinates(Ve,1);
        ye = geom.elements.coordinates(Ve,2);
        
        ii = geom.pivot.pivot(Vb);
        if ii > 0 % se fosse negativo sarebbe di Dirichlet (ha priorità)  
            temp = 0;
            for k = 1:5
                temp = temp + ni*w(k)*gn_fun(xb+(xe-xb)*t_hat(k), yb+(ye-yb)*t_hat(k), t, geom.pivot.Ne(b,2))*phib(k)*E_norm;
            end
            gn(ii) = gn(ii) + temp; 
        end

        ii = geom.pivot.pivot(Ve);
        if ii > 0   
            temp = 0;
            for k = 1:5
                temp = temp + ni*w(k)*gn_fun(xb+(xe-xb)*t_hat(k), yb+(ye-yb)*t_hat(k), t, geom.pivot.Ne(b,2))*phie(k)*E_norm;
            end
            gn(ii) = gn(ii) + temp; 
        end


        ii = pivotMid(pivotMidNe(b,1),1);  % prima dato un lato avevi due nodi, ora dato un lato stai già identificando quale nodo interno stai considerando
        temp=0;
        for k = 1:5
            temp = temp + ni*w(k)*gn_fun(xb+(xe-xb)*t_hat(k), yb+(ye-yb)*t_hat(k), t, pivotMidNe(b,2))*phim(k)*E_norm;
        end
        gn(ii+sum(geom.pivot.pivot >0)) = gn(ii+sum(geom.pivot.pivot >0)) + temp;       
    end
    

    

end