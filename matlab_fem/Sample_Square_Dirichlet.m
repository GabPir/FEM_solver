    
function [] =  Sample_Square_Dirichlet(area_max, indice_figura, stampa)

    global geom;


    Domain.InputVertex = [ 0 0
                       1 0
                       1 1
                       0 1];


    Domain.Boundary.Values = 1:4;   % definiamo 4 lati che prende automaticamente
    Domain.Holes.Hole = [];       % non ci sono buchi nel dominio
    Domain.Segments.Segment = []; % non ci sono lati forzati nel dominio
    
    % valori numerici per le condizioni al contorno 
    BC.Values = [0.0 12.0 0.0 14.0 0.0 16.0 0.0 0.0 0.0];
    
    %Marker dispari -> Dirichlet; pari -> Neumann (dovrà metterci gradi di libertà data l'assenza di valori della funzioni 
    BC.Boundary.Values = [2 1 3 5];  
    % marker dei Vertici iniziali, ad ogni nodo associo un marker dispari generico che sia DIVERSO DAGLI ALTRI NODI (agli alti nodi darà il valore che si ha associato al lato dove risiede il nodo (3,5,7,9)
    BC.InputVertexValues = [5 1 3 5];
    

    
    
    BC.Holes.Hole = [];   % non ci sono stranezze nelle condizioni (omogenee)
    BC.Segments.Segment = [];


    
    
    % Inserimento dei parametri di triangolazione
    RefiningOptions.CheckArea  = 'Y'; %chiedo  di costgruire triangoli cui aria è inferiore a quella che gli prescrivo
    RefiningOptions.CheckAngle = 'N'; %assicuro che l'angolo minimo sia superiore ad un certo valore (mesh di buona qualità) -> nel quadrato non fa differenza
    RefiningOptions.AreaValue  = area_max; % se avessi messo YES sopra metto un ancolo minimo (mesh di qualità è meglio avere come minmimo 20)
    RefiningOptions.AngleValue = []; %prescrive un meshamento diverso per diverse regioni
    RefiningOptions.Subregions = [];
    
    [geom] = bbtr30(Domain,BC,RefiningOptions);   % TRIANGOLATORE: disegna la griglia 
    if stampa == 1
        draw_grid (geom,indice_figura); 
    end
    
    %Allocazione corretta memoria e pivot (gradi di liberà)
    
    geom.elements.coordinates = geom.elements.coordinates(...
				    1:geom.nelements.nVertexes,:);
    geom.elements.triangles = geom.elements.triangles(...
				    1:geom.nelements.nTriangles,:);
    geom.elements.borders = geom.elements.borders(...
				    1:geom.nelements.nBorders,:);
    geom.elements.neighbourhood = geom.elements.neighbourhood(...
				    1:geom.nelements.nTriangles,:);
    
    
    
    
    j  = 1;
    Dj = 1;
    for i=1:size(geom.pivot.nodelist)
         if geom.pivot.nodelist(i)==0
            geom.pivot.pivot(i)=j;
            j = j+1;
         else
            geom.pivot.pivot(i)=-Dj;
            Dj = Dj + 1;
         end
    end
    
    geom.pivot.pivot = transpose(geom.pivot.pivot);
    
    [X,I] = sort(geom.pivot.Di(:,1));
    geom.pivot.Di = geom.pivot.Di(I,:);
    
    clear X I;
    






end
