function [] =  Numdof() 
   global geom Ndof NdMid Nd;

   Ntot = sum(geom.pivot.pivot >0);
   for i = 1:geom.nelements.nBorders
        if geom.elements.borders(i,3) == -1 || geom.elements.borders(i,4)==-1
            if geom.pivot.pivot(geom.elements.borders(i,1)) > 0 || geom.pivot.pivot(geom.elements.borders(i,2)) > 0
                Ntot = Ntot +1;
            end
       
        else
            Ntot = Ntot +1;
        end
   end

   NdMid = 0;
   for i = 1:geom.nelements.nBorders
        if geom.elements.borders(i,3) == -1 || geom.elements.borders(i,4)==-1
            if geom.pivot.pivot(geom.elements.borders(i,1)) < 0 && geom.pivot.pivot(geom.elements.borders(i,2)) < 0
                NdMid = NdMid +1;
            end
        end
   end


   Ndof = Ntot;
   NdMid = NdMid;
   Nd = sum(geom.pivot.pivot < 0) + NdMid; 
