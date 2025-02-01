function [f_0, gd_0, gd_0_prime, x_0, gn_0] = CondizioneIniziale(SUPG)
    global u0 geom PivotMid MidCoordinates pivotMid pivotMidDi NdMid Ndof Nd t gd gd_prime gn f forzante;
    
    f_0 = zeros(Ndof, 1);   
    gd_0 = zeros(Nd,1);
    gd_0_prime = zeros(Nd,1);
    x_0 = zeros(Ndof, 1);
    gn_0 = zeros(Nd,1);
    
    for i=1:geom.nelements.nVertexes
        ii = geom.pivot.pivot(i);
        if ii > 0
            x_0(ii) = uIniziale(geom.elements.coordinates(i,1), geom.elements.coordinates(i,2), t); % HAI COSTRUITO PIVOT IN MANIERA TALE CHE LE POSIZIONI CORRISPONDANO AGLI INDICI GLOBALI
        end
    end

    for i=1:length(pivotMid)
        ii = pivotMid(i);
        if ii > 0
            x_0(ii+sum(geom.pivot.pivot()>0)) = uIniziale(MidCoordinates(i,1), MidCoordinates(i,2), t);  % la posizione in x non è shiftata di # vertici ma solo dei Ndof di questi ultimi
        end
    end


    % la condizione iniziale deve soddisfare tutte le condizioni a t = 0:
    AggiungiForzante(SUPG);
    AggiungiNeumann();
    AggiungiDirichlet();
    AggiungiDirichletDerivata();
    gd_0 = gd;
    gn_0 = gn;
    gd_0_prime = gd_prime;
    f_0 = f;
end

