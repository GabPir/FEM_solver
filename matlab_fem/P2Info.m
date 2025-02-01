function [] = P2Info()
    global geom NdMid;
    global EleMid pivotMid pivotMidDi pivotMidNe MidCoordinates;
    EleMid = zeros(geom.nelements.nTriangles, 3);
    pivotMid = zeros(geom.nelements.nBorders, 1);
    pivotMidDi = zeros(NdMid,1);
    pivotMidNe = zeros(length(geom.pivot.Ne(:,1)), 1);
    MidCoordinates = zeros(geom.nelements.nBorders, 2);
    j = -1;
    k = 1;
    t = 1;
    for i = 1:geom.nelements.nBorders
        Vbe = geom.elements.borders(i,1);
        Vfi = geom.elements.borders(i,2);
        T1 = geom.elements.borders(i,3);
        T2 = geom.elements.borders(i,4);

        MidCoordinates(i, 1) = geom.elements.coordinates(Vbe,1) + 0.5*(geom.elements.coordinates(Vfi,1)-geom.elements.coordinates(Vbe,1));
        MidCoordinates(i, 2) = geom.elements.coordinates(Vbe,2) + 0.5*(geom.elements.coordinates(Vfi,2)-geom.elements.coordinates(Vbe,2));


   
        if T1 > 0
            v = geom.elements.triangles(T1,:);
            if isequal([v(2),v(3)],[Vbe,Vfi]) || isequal([v(2),v(3)],[Vfi,Vbe])
                EleMid(T1,3) = i;
            elseif isequal([v(2),v(1)], [Vbe,Vfi]) || isequal([v(2),v(1)],[Vfi,Vbe])
                EleMid(T1,2) = i;
            else
                EleMid(T1,1) = i;
            end
        end 
        if T2 > 0
            v = geom.elements.triangles(T2,:);
            if isequal([v(2),v(3)],[Vbe,Vfi]) || isequal([v(2),v(3)],[Vfi,Vbe])
                EleMid(T2,3) = i;
            elseif isequal([v(2),v(1)], [Vbe,Vfi]) || isequal([v(2),v(1)],[Vfi,Vbe])
                EleMid(T2,2) = i;
            else
                EleMid(T2,1) = i;
            end
        end
        
        if rem(geom.support.BInfo(i,3),2) == 0 || geom.support.BInfo(i,3) == 0
            pivotMid(i,1) = k;
            if geom.support.BInfo(i,3) > 0
                pivotMidNe(t,1) = i;
                pivotMidNe(t,2) = geom.support.BInfo(i,3);
                t = t +1;
            end
            k = k+1;
        elseif rem(geom.support.BInfo(i,3),2) ~= 0
            pivotMid(i,1) = j;
            pivotMidDi(-j,1) = i;
            pivotMidDi(-j,2) = geom.support.BInfo(i,3);
            j = j -1;
        end

    end 
end