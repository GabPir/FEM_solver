function [value] = QuadraturaErroreH1_P2(e, grad_u_esatta)
    global geom EleMid t;
    global u;
    
    % valori formula di quadratura
    [x_hat,y_hat,~,w] = int_nodes_weights(5);
    x_hat=x_hat';
    y_hat=y_hat';


    x = geom.elements.coordinates(:,1);
    y = geom.elements.coordinates(:,2);
    
   
    a1 = geom.elements.coordinates(geom.elements.triangles(e,1),:)';
    a2 = geom.elements.coordinates(geom.elements.triangles(e,2),:)';
    a3 = geom.elements.coordinates(geom.elements.triangles(e,3),:)';

    Ele = [geom.elements.triangles, (EleMid+geom.nelements.nVertexes)];

    
    b11 = x(Ele(e,1))-x(Ele(e,3)); %x(Ele(e,2))-x(Ele(e,1));
    b12 = x(Ele(e,2))-x(Ele(e,3)); %x(Ele(e,3))-x(Ele(e,1));
    b21 = y(Ele(e,1))-y(Ele(e,3)); %y(Ele(e,2))-y(Ele(e,1));
    b22 = y(Ele(e,2))-y(Ele(e,3)); %y(Ele(e,3))-y(Ele(e,1)); 
    B = [b11 b12; b21 b22];
    BinvT = inv(B');


    GradPhi = [];  %nodo/derivata parziale/base
    GradPhi(:,:,1) = [4*(x_hat(:,1))-1, zeros(7,1)];
    GradPhi(:,:,2) = [zeros(7,1), 4*(y_hat(:,1))-1];
    GradPhi(:,:,3) = [-3+4*x_hat(:,1)+4*y_hat(:,1), -3+4*x_hat(:,1)+4*y_hat(:,1)];
    GradPhi(:,:,4) = [-8*x_hat(:,1)-4*y_hat(:,1)+4,-4*x_hat(:,1)];
    GradPhi(:,:,5) = [4*y_hat(:,1),4*x_hat(:,1)];
    GradPhi(:,:,6) = [-4*y_hat(:,1),-8*y_hat(:,1)-4*x_hat(:,1)+4];
    

    value = 0;
    for k = 1:7 %Nq
        coord = a3.*(1-x_hat(k)-y_hat(k)) + a1.*x_hat(k) + a2.*y_hat(k); % formula del cambio variabili
        vet_temp = zeros(2,1);
        for i = 1:6 
            vet_temp = vet_temp + u(Ele(e,i))*GradPhi(k,:,i)';
        end
        value = value + w(k)*(grad_u_esatta([coord',t])' - BinvT*vet_temp)'*(grad_u_esatta([coord',t])' - BinvT*vet_temp);
    end 
    
