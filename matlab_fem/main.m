 global geom ni sigma Beta;
 global A A_d f u gd gd_prime gn Nd NdMid Ndof C C_d R R_d t M M_d K K_d;
 global EleMid pivotMid pivotMidDi pivotMidNe MidCoordinates forzante;


 area_min = [0.005];
 Ns = [10 50 100];
 area_min = flip(area_min);

 Beta = [1,0];
 sigma = 0.01;
 ni = 0.01;

 SUPG = 0; % 1 = attivo 0 = inattivo
 Mass_lumping = 0;
 
 forzante = @(v) cos(v(:,3)) -ni*4;
 u_esatta = @(v) sin(v(:,3)) + v(:,1).^2 + v(:,2).^2;
 grad_u_esatta = @(v) Grad_u_esatta(v);
 

 eL2 = zeros(length(area_min),length(Ns));
 eLinfinity = zeros(length(area_min),length(Ns));
 eH1 = zeros((length(area_min)),length(Ns));

indice_figura = 0;

for j = 1:length(Ns)
     N = Ns(j);
     dt = 1/N;
     times = linspace(0,1,N+1);
     indice_figura = 0;
     for k = 1:length(area_min)
        
        indice_figura = indice_figura + 1;
        Sample_Square_Dirichlet(area_min(k), indice_figura, 1); % 1 se vuoi stampare 0 altrimenti
        Numdof();
        P2Info();
        
    
        A = spalloc(Ndof, Ndof, 10*Ndof);
        C = spalloc(Ndof, Ndof, 10*Ndof);
        R = spalloc(Ndof, Ndof, 10*Ndof);
        A_d = spalloc(Ndof, Nd, 10*Ndof); 
        C_d = spalloc(Ndof, Nd, 10*Ndof);
        R_d = spalloc(Ndof, Nd, 10*Ndof);
        f = zeros(Ndof, 1);
        u = zeros(Ndof + Nd, 1);
        gd = zeros(Nd,1);
        gd_prime = zeros(Nd,1);
        gn = zeros(Ndof, 1);
    
        GalerkinDC(SUPG,Mass_lumping)   % SUPG attivo 1, Mass Lumping attivo 1 
        
        t=0;
        [f_old, gd_old, gd_prime_old, x_old, gn_old] = CondizioneIniziale(SUPG);
        AssemblaSoluzione(x_old)  
        StampaSoluzione(indice_figura);
        
        for index = 2:length(times)
            t = times(index);
            f = zeros(Ndof, 1);
            u = zeros(Ndof + Nd, 1);
            gd = zeros(Nd,1);
            gd_prime = zeros(Nd,1);
            gn = zeros(Ndof, 1);
            
            AggiungiForzante(SUPG)
            AggiungiNeumann()
            AggiungiDirichlet()
            AggiungiDirichletDerivata()
            
            x = (M+(dt/2)*K)\((M-(dt/2)*K)*x_old-(dt/2)*(M_d*gd_prime_old + K_d*gd_old + M_d*gd_prime +K_d*gd -(f_old+f+gn_old+gn))); 
            AssemblaSoluzione(x)
            
            
            if index == round(length(times)/2) || index == length(times)
                indice_figura = indice_figura + 1;
                StampaSoluzione(indice_figura);
            end
    
            if index == length(times)
                [eLinfinity(k,j), eL2(k,j), eH1(k,j)]  = StimaErroreP2(u_esatta, grad_u_esatta);
            end
    
            gn_old = gn;
            gd_old = gd;
            gd_prime_old = gd_prime;
            f_old = f;
            x_old = x;
        end 
    end
end 



log_eLinfinity = log(eLinfinity(end,:)');
log_eL2 = log(eL2(end,:)');
log_eH1 = log(eH1(end,:)');
log_area_min = log(area_min);
indice_figura = indice_figura+1;
%figure(indice_figura)
%plot(log_area_min, log_eLinfinity)
indice_figura = indice_figura+1;
figure(indice_figura)
plot(log([1/Ns(1); 1/Ns(2); 1/Ns(3)]), log_eL2)
title("L2")
indice_figura = indice_figura+1;
figure(indice_figura)
plot(log([1/Ns(1); 1/Ns(2); 1/Ns(3)]),log_eH1)
title("H1")

disp("valori fitting logaritmico errore L2:")
polyfit(log([1/Ns(1); 1/Ns(2); 1/Ns(3)]), log_eL2,1)
disp("valori fitting logaritmico errore H1:")
polyfit(log([1/Ns(1); 1/Ns(2); 1/Ns(3)]), log_eH1,1) 


% disp("valori fitting logaritmico errore L2:")
% polyfit(log((sqrt(area_min))'), log_eL2,1)
% disp("valori fitting logaritmico errore H1:")
% polyfit(log((sqrt(area_min))'), log_eH1,1)  



function [d_v] = Grad_u_esatta(v)
    d_v(:,1) = 2*v(:,1);
    d_v(:,2) = 2*v(:,2);
end
 