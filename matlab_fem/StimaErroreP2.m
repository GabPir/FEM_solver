function [eLinfty, eL2, eH1] = StimaErroreP2(u_esatta, grad_u_esatta)
    global geom;
    global u MidCoordinates t;
    
    % Stima errore in norma L  infinito:
    XY = [geom.elements.coordinates;MidCoordinates];
    XYT = [XY, t*ones(length(XY),1)];
    eLinfty = max(abs(u_esatta(XYT) - u));
    
    
    % Stima errore in norma L 2:    
    eL2 = 0;
    for e = 1:geom.nelements.nTriangles
        A = geom.support.TInfo(e).Area;
        eL2 = eL2 + 2*A*QuadraturaErroreL2_P2(e, u_esatta);
    end 
    eL2 = (eL2)^0.5;

    % Stima errore in H1: 
    % E' EQUIVALENTE USARE SOLO LA COMPONENTE DEL VETTORE
    eH1 = 0;
    for e = 1:geom.nelements.nTriangles
        A = geom.support.TInfo(e).Area;
        eH1 = eH1 + 2*A*QuadraturaErroreH1_P2(e, grad_u_esatta); 
    end 
    eH1 = (eH1)^0.5;
end