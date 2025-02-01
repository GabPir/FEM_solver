function [] = GalerkinDC(SUPG, Mass_lumping)

    global geom ni Beta sigma;
    global A A_d gd C C_d R R_d Ndof Nd  M_d M K K_d;
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
    

    GradPhi_xx = [];
    GradPhi_xx(1,:) = 4*ones(1,7);
    GradPhi_xx(2,:) = 0*ones(1,7);
    GradPhi_xx(3,:) = 4*ones(1,7);
    GradPhi_xx(4,:) = -8*ones(1,7); 
    GradPhi_xx(5,:) = 0*ones(1,7); 
    GradPhi_xx(6,:) = 0*ones(1,7);


    GradPhi_xy = [];
    GradPhi_xy(1,:) = 0*ones(1,7);
    GradPhi_xy(2,:) = 0*ones(1,7);
    GradPhi_xy(3,:) = 4*ones(1,7);
    GradPhi_xy(4,:) = -4*ones(1,7); 
    GradPhi_xy(5,:) = 4*ones(1,7); 
    GradPhi_xy(6,:) = -4*ones(1,7);

    GradPhi_yy = [];
    GradPhi_yy(1,:) = 0*ones(1,7);
    GradPhi_yy(2,:) = 4*ones(1,7);
    GradPhi_yy(3,:) = 4*ones(1,7);
    GradPhi_yy(4,:) = 0*ones(1,7); 
    GradPhi_yy(5,:) = 0*ones(1,7); 
    GradPhi_yy(6,:) = -8*ones(1,7);



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
            if ii > 0  % se fosse negativo sarebbe un nodo al bordodi dirichlet (no gradi di libertà)
                if i > 3
                    ii = ii +max(geom.pivot.pivot);  % Riempi la matrice A da dopo le posizioni nei nodi
                end

                for j=1:6
                    if j <= 3
                        jj = pivot(Ele(e,j));
                    else
                        jj = pivotMid(Ele(e,j));
                    end
                    if jj > 0
                        if j > 3
                            jj = jj + max(geom.pivot.pivot); % Riempi la matrice A da dopo le posizioni nei nodi
                        end
                        
                        temp = 0;
                        temp2 = 0;
                        temp3 = 0;
                        for k = 1:7
                            temp = temp + ni*Area*2*GradPhi(k,:,j)*T*GradPhi(k,:,i)'*w(k);
                            temp2 = temp2 + Beta*inv(B')*GradPhi(k,:,j)'*Phi(k,i)*2*Area*w(k);
                            temp3 = temp3 + sigma*Phi(k,i)*Phi(k,j)*2*Area*w(k);
                            if Peclet > 1
                                temp = temp - Tau*ni*Area*2*(T(1,:)*[GradPhi_xx(j,k); GradPhi_xy(j,k)] + T(2,:)*[GradPhi_xy(j,k); GradPhi_yy(j,k)])*Beta*inv(B')*GradPhi(k,:,i)'*w(k);
                                temp2 = temp2 + Tau*Area*2*Beta*inv(B')*GradPhi(k,:,j)'*Beta*inv(B')*GradPhi(k,:,i)'*w(k);
                            end
                        end
                        A(ii,jj) = A(ii,jj) + temp;
                        C(ii,jj) = C(ii,jj) + temp2;
                        R(ii,jj) = R(ii,jj) + temp3;
                    else
                        % DIRICHLET NON OMOGENEO:
                        if j > 3
                            jj = jj - sum(geom.pivot.pivot(:,1)<0);  
                        end
                        temp = 0;
                        temp2 = 0;
                        temp3 = 0;
                        for k = 1:7
                            temp = temp + ni*Area*2*GradPhi(k,:,j)*T*GradPhi(k,:,i)'*w(k);
                            temp2 = temp2 + Beta*inv(B')*GradPhi(k,:,j)'*Phi(k,i)*2*Area*w(k);
                            temp3 = temp3 + sigma*Phi(k,i)*Phi(k,j)*2*Area*w(k);
                            if Peclet > 1 
                                temp = temp - Tau*ni*Area*2*(T(1,:)*[GradPhi_xx(j,k); GradPhi_xy(j,k)] + T(2,:)*[GradPhi_xy(j,k); GradPhi_yy(j,k)])*Beta*inv(B')*GradPhi(k,:,i)'*w(k); 
                                temp2 = temp2 + Tau*Area*2*Beta*inv(B')*GradPhi(k,:,j)'*Beta*inv(B')*GradPhi(k,:,i)'*w(k);
                            end
                        end
                        A_d(ii,-jj) = A_d(ii,-jj) + temp; 
                        C_d(ii,-jj) = C_d(ii,-jj) + temp2;
                        R_d(ii,-jj) = R_d(ii,-jj) + temp3;
                        
                    end 

                end

             end 
     
         end
    end
    
    if Mass_lumping == 1
        Temp = zeros(Ndof,Ndof);
        for i = 1:Ndof
            Temp(i,i)=sum(R(i,:));
        end 
        R = Temp;
    end

    

    K = A;
    K_d = A_d;
    M = (1/sigma)*R; %la matrice relativa al termine temporale si construisce come R
    M_d = (1/sigma)*R_d;
    
    