function [] =  AssemblaSoluzione(x)
    global geom u pivotMid MidCoordinates pivotMidDi t;

    for i=1:geom.nelements.nVertexes
        ii = geom.pivot.pivot(i);
        if ii > 0
            u(i) = x(ii); % HAI COSTRUITO PIVOT IN MANIERA TALE CHE LE POSIZIONI CORRISPONDANO AGLI INDICI GLOBALI
        else
            u(i) = gd_fun(geom.elements.coordinates(i,1), geom.elements.coordinates(i,2),t, geom.pivot.Di(-ii,2)); % condizione di dirichlet al contorno 
        end
    end

    for i=1:length(pivotMid)
        ii = pivotMid(i);
        if ii > 0
            u(i+geom.nelements.nVertexes) = x(ii+sum(geom.pivot.pivot > 0));  % la posizione in x non è shiftata di # vertici ma solo dei Ndof di questi ultimi
        else
            u(i+geom.nelements.nVertexes) = gd_fun(MidCoordinates(i,1), MidCoordinates(i,2),t, pivotMidDi(-ii,2)); 
        end
    end
end