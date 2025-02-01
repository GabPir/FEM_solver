function [] = AggiungiDirichletDerivata()
    global geom gd_prime Nd NdMid MidCoordinates t pivotMidDi;
    
    for i=1:Nd-NdMid
        x = geom.elements.coordinates(geom.pivot.Di(i,1),1);
        y = geom.elements.coordinates(geom.pivot.Di(i,1),2);
        gd_prime(i) = gd_fun_prime(x, y, t,  geom.pivot.Di(i,2)); 

    end          
    xmid = MidCoordinates(:,1);
    ymid = MidCoordinates(:,2);
    for j=1:NdMid
        k = (Nd-NdMid)+j;
        gd_prime(k) = gd_fun_prime(xmid(pivotMidDi(j,1)), ymid(pivotMidDi(j,1)), t, pivotMidDi(j,2));
    end
end