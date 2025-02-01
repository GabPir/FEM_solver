function [value] = QuadraturaErroreL2_P2(e, u_esatta, t)
    global geom EleMid t;
    global u MidCoordinates;
    
    % valori formula di quadratura
    [x_hat,y_hat,~,w] = int_nodes_weights(5);


    a1 = geom.elements.coordinates(geom.elements.triangles(e,1),:)';
    a2 = geom.elements.coordinates(geom.elements.triangles(e,2),:)';
    a3 = geom.elements.coordinates(geom.elements.triangles(e,3),:)';
    value = 0;


    for k = 1:6 %Nq
        % Le phy valutate sui punti sono direttamente scritte in formula 
        coord = a3.*(1-x_hat(k)-y_hat(k)) + a1.*x_hat(k) + a2.*y_hat(k); % formula del cambio variabili
        temp= u_esatta([coord',t])...
               -(u(geom.elements.triangles(e,1))*2*x_hat(k)*(x_hat(k)-0.5)...
                + u(geom.elements.triangles(e,2))*2*y_hat(k)*(y_hat(k)-0.5)...
                + u(geom.elements.triangles(e,3))*2*(1-x_hat(k)-y_hat(k))*(0.5-x_hat(k)-y_hat(k))...
                + u(geom.nelements.nVertexes+EleMid(e,1))*4*(1-x_hat(k)-y_hat(k))*x_hat(k)...
                + u(geom.nelements.nVertexes+EleMid(e,2))*4*x_hat(k)*y_hat(k)...
                + u(geom.nelements.nVertexes+EleMid(e,3))*4*(1-x_hat(k)-y_hat(k))*y_hat(k)); 
        
        value = value + w(k)*(temp)^2;
    end 
    
