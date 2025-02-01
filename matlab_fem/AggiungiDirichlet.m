function [] = AggiungiDirichlet()
    global geom gd Nd NdMid MidCoordinates pivotMidDi t;
    
    for i=1:Nd-NdMid
        x = geom.elements.coordinates(geom.pivot.Di(i,1),1);
        y = geom.elements.coordinates(geom.pivot.Di(i,1),2);
        gd(i) = gd_fun(x, y, t,  geom.pivot.Di(i,2)); 

    end          
    xmid = MidCoordinates(:,1);
    ymid = MidCoordinates(:,2);
    for j=1:NdMid
        k = (Nd-NdMid)+j;
        
        gd(k) = gd_fun(xmid(pivotMidDi(j,1)), ymid(pivotMidDi(j,1)), t,  pivotMidDi(j,2));
    end
end