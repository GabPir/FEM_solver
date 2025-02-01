function StampaSoluzione(indice_figura)
    global geom EleMid MidCoordinates u;
    figure(indice_figura)
    T = [geom.elements.triangles; EleMid+geom.nelements.nVertexes];
    x = [geom.elements.coordinates(:,1); MidCoordinates(:,1)];
    y = [geom.elements.coordinates(:,2); MidCoordinates(:,2)];
    trisurf(T, x, y, u)
end 