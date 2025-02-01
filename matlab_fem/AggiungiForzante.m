function [] = AggiungiForzante(SUPG)

    global geom ni Beta sigma ni forzante;
    global A A_d f gd C C_d t R R_d Ndof Nd;
    global EleMid pivotMid pivotMidDi pivotMidNe MidCoordinates;
    
    
    
        
    [x_hat,y_hat,~,w] = int_nodes_weights(5);
    x_hat=x_hat';
    y_hat=y_hat';
    Phi = [2.*x_hat(:,1).*(x_hat(:,1)-1/2), ...
           2.*y_hat(:,1).*(y_hat(:,1)-1/2), ... 
           2.*(1-x_hat(:,1)-y_hat(:,1)).*(1-x_hat(:,1)-y_hat(:,1)-1/2), ...
           4.*(1-x_hat(:,1)-y_hat(:,1)).*x_hat(:,1), ...
           4.*y_hat(:,1).*x_hat(:,1), ...
           4.*y_hat(:,1).*(1-x_hat(:,1)-y_hat(:,1))];
    
    
    GradPhi = [];  %nodo/derivata parziale/base
    GradPhi(:,:,1) = [4*(x_hat(:,1))-1, zeros(7,1)];
    GradPhi(:,:,2) = [zeros(7,1), 4*(y_hat(:,1))-1];
    GradPhi(:,:,3) = [-3+4*x_hat(:,1)+4*y_hat(:,1), -3+4*x_hat(:,1)+4*y_hat(:,1)];
    GradPhi(:,:,4) = [-8*x_hat(:,1)-4*y_hat(:,1)+4,-4*x_hat(:,1)];
    GradPhi(:,:,5) = [4*y_hat(:,1),4*x_hat(:,1)];
    GradPhi(:,:,6) = [-4*y_hat(:,1),-8*y_hat(:,1)-4*x_hat(:,1)+4];
    
    Ele = [geom.elements.triangles, EleMid]; %matrice T con anche i nodi interni 4, 5, 6
    Nt = geom.nelements.nTriangles; %numero triangoli
    x = geom.elements.coordinates(:,1);
    y = geom.elements.coordinates(:,2);
    pivot = geom.pivot.pivot;

    
    for e=1:Nt
        %calcoliamo i valori geometrici
        
        Area = geom.support.TInfo(e).Area;
        b11 = x(Ele(e,1))-x(Ele(e,3)); %x(Ele(e,2))-x(Ele(e,1));
        b12 = x(Ele(e,2))-x(Ele(e,3)); %x(Ele(e,3))-x(Ele(e,1));
        b21 = y(Ele(e,1))-y(Ele(e,3)); %y(Ele(e,2))-y(Ele(e,1));
        b22 = y(Ele(e,2))-y(Ele(e,3)); %y(Ele(e,3))-y(Ele(e,1)); 
        B = [b11 b12; b21 b22];
        T = inv(B)*inv(B');
         
        h = max([norm(geom.elements.coordinates(Ele(e,1),:)-geom.elements.coordinates(Ele(e,2),:),2), ...
            norm(geom.elements.coordinates(Ele(e,2),:)-geom.elements.coordinates(Ele(e,3),:),2), ...
            norm(geom.elements.coordinates(Ele(e,3),:)-geom.elements.coordinates(Ele(e,1),:),2)]);
                
        Peclet = (1/(2*24*ni))*norm(Beta,2)*h;
        Tau = 0;
        if SUPG == 1
            if Peclet > 1
                Tau = h/(2*norm(Beta,2));
            else
                Tau = (h^2)/(24*4*ni);
            end
        end

        for i=1:6
            if i <= 3
                ii = pivot(Ele(e,i));
            else
                ii = pivotMid(Ele(e,i));
            end
            if ii > 0 
                if i > 3
                    ii = ii +max(geom.pivot.pivot);  
                end
                temp = 0;
                for k = 1:7
                    vett = [x(Ele(e,3)), y(Ele(e,3))]'+B*[x_hat(k), y_hat(k)]';
                    temp = temp + 2*Area*w(k)*forzante([vett(1), vett(2), t])*Phi(k,i);
                    if Peclet > 1
                        temp = temp + Tau*2*Area*w(k)*forzante([vett(1), vett(2), t])*Beta*inv(B')*GradPhi(k,:,i)';
                    end
                end
                f(ii) = f(ii) + temp;
             end 
         end
    
    end