//Implementation of a clock model (neurospora-mice)

//Fonction de Hill

function y = hillp(x, K, n)
    y = x^n / (K^n + x^n);
endfunction


function v=clock_model_rhs(t,x,par)
    
//Variables
    W=x(1); //Whitecollar, activator (Bmal1)
    frq=x(2); //frq mRNA (per)
    Pcy=x(3); //FRQ prot, cytoplasme (PER)
    Pnu=x(4); //FRQ prot, nucleus
    
//Parameters
    a1=par(1);
    a2=par(2);
    a3=par(3);
    a4=par(4);
    gamma=par(5);
    gM=par(6);
    gP=par(7);
    
    K1=par(8);
    K4=par(9);
    
    n=par(10);
    
    deg=gamma*Pnu*W;
    
    v=[];
    
    v(1)=a4*hillp(W,K4,n)-deg-gP*W;
    v(2)=a1*hillp(W,K1,n)-gM*frq;
    v(3)=a2*frq-gP*Pcy-a3*Pcy;
    v(4)=a3*Pcy-gP*Pnu-deg;
endfunction

//Test avec des paramètres dans les intervalles
a1 = 8;
a2 = 7;
a3 = 0.7;
a4 = 40;
gamma = 0.15;
gP = 0.01; // Correspond à gamma P
gM = 0.4; // doit être inférieur à  1/2 * a4/K4
K1 = 50;
K4 = 0.7;
n = 2;

par = [a1, a2, a3, a4, gamma, gM, gP, K1, K4, n];

x0 = [0.1; 1; 1; 5]; // Conditions initiales
t0 = 0;
dt = 0.01;
tvec = t0:dt:150;





//Recherche des points d'équilibre

//Graphiquement

W=[0:0.01:100];
for i=1:length(W)
    g1(i)=a4*hillp(W(i),K4,n)/(gamma*W(i)) -gP/gamma;
    aux=a3/(gP+gamma*W(i));
    g2(i)=(a2/(a3+gP))*aux*(a1/gM)*hillp(W(i),K1,n);
end
f3=figure(3);
plot(W,g1,"r-",W,g2,"b-",'linedwidth',3);
f3.background=color("white");
legend('dW/dt=0','dPnu/dt=0');



//Numériquement

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
function Weq = Wintersection(Wo, par)
// cette fonction trouve les zeros associés à l'équation 
// f= g1(0)-g2(0)
    a1 = par(1);
    a2 = par(2);
    a3 = par(3);
    a4 = par(4);
    gamma = par(5);
    gM = par(6);
    gP = par(7); 
    n = par(8);
    K1 = par(9); 
    K4 = par(10); 
    aux1 = a3 / (gamma * Wo + gP);
    aux2 = a2 / (a3 + gP);
    aux3 = (a1 / gM) * ((Wo^n) / (K1^n + Wo^n)); 
    g1 = aux1 * aux2 * aux3;
    g2 = (a4 * (Wo^(n-1)) / (K4^n + Wo^n) - gP) / gamma; 
    Weq = g2-g1 ; 
endfunction


//Definir W0 'initial guest'
W0 = 4 * K4; // Supposant que vous voulez multiplier K4 par 4

[Wst, feq] = fsolve(Wintersection, W0, par);

// Définition de W0 ou "initial guess" Et définition de tous les points d'équilibre 
//W0 = K4 * 4;
//[W_st, feq] = fsolve(W0, list(W_intersection, para)); // K4 à la place de W0 ? K4 est proche de W0
//m_st = (a1/gM) / (W^n + K1^n);
//Pc_st = (a2 / (a3 + gP)) * m_st;                     JSP SI IL FAUT METTRE CA 
//Pn_st = (a3 / (gP + gamma * W_st)) * Pc_st;

// Matrice jacobienne

A = [];
// dW/dt
A(1,1) = (a4 * n * W_st^(n-1) / (K4^n + W_st^n)^2) -gP -gamma * Pn_st; // Dérivée partielle des équations de base par rapport à W
A(1,2) = 0; // Par rapport à m
A(1,3) = 0; // Par rapport à Pc
A(1,4) = -gamma * W_st; // Par rapport à Pn
// 2e ligne matrice / 2e équation dm/dt
A(2,1) = (a1 * n * K1^n * W_st^(n-1)) / (K1^n + W_st^n)^2;
A(2,2) = -gM;
A(2,3) = 0;
A(2,4) = 0;
// dPc/dt
A(3,1) = 0;
A(3,2) = a2;
A(3,3) = - (a3 + gP);
A(3,4) = 0;
// dPn/dt
A(4,1) = -gamma * Pn_st;
A(4,2) = 0;
A(4,3) = a3;
A(4,4) = -gP -gamma * W_st;

// Calcul des valeurs propres partie réelle négative --> poit équilibre stables
vp = spec(A);

// Calcul matrice jacobienne à W_st = Pc_st = Pn_st = m_st = 0
A0 = [];

A0(1,1) = -gP;
A0(1,2) = 0; 
A0(1,3) = 0; 
A0(1,4) = 0;
// 2e ligne matrice / 2e équation dm/dt
A0(2,1) = 0;
A0(2,2) = -gM;
A0(2,3) = 0;
A0(2,4) = 0;
// dPc/dt
A0(3,1) = 0;
A0(3,2) = a2;
A0(3,3) = - (a3 + gP);
A0(3,4) = 0;
// dPn/dt
A0(4,1) = 0;
A0(4,2) = 0;
A0(4,3) = a3;
A0(4,4) = -gP;

//Pour calculer la période

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
function Periode=calculperiode(X,tvec)
//Initialisation
Periode = 0 ;
j=0;
//X est un vecteur X=[X1,...,Xn]
L=length(X);
S=floor(0.3*L);
MaxX=max(X(S:L));
MinX=min(X(S:L));
theta=0.5*(MaxX+MinX)
for i=S : L-1 
    auxg=X(i)-theta;
    auxd=X(i+1)-theta; 
    if and ([auxg<0 , auxd>0])
        j=j+1;
        Pvec(j)=tvec(i);
    end 
end

//Pvec contient les instants de tous les croisements(t1,t2,t3,etc)
Periode=mean(diff(Pvec));
//diff(Pvec)= [Pvec(2)-Pvec(1),Pvec(3)-Pvec(2),....,Pvec(n)-Pvec(n-1)] 
endfunction

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
function S=sensitivity(par, tvec)
    //Calcul de la période "original"
    t0=0, X0=[1;5;5;1];
    
    sol=ode("stiff",x0,t0,tvec, list(clock_model_rhs,par) );
    Periode=calculperiode(sol(1,:),tvec);
    
    //Calcul des sensibilités
    for j=1:10
        parS=par;
        parS(j)=parS(j)+0.05*parS(j);
        solS=ode("stiff",x0,t0,tvec, list(clock_model_rhs,parS) );
        PerS=calculperiode(solS(1,:),tvec);
        S(j)=(PerS-Periode)/(parS(j)-par(j));
    end
    Smean=mean(S.*S);
endfunction



//Partie principale

sol=ode("stiff",x0,t0,tvec, list(clock_model_rhs,par) );
Periode=calculperiode(sol(1,:),tvec),
//Faire le graphe de la solution au cour du temps, dans la figure 1
//Options: couleur (r=red, b=blue,k=black,...), forme de ligne (-,--,:,-.)
figure(1);
plot(tvec,sol(1,:),'b-',tvec,sol(2,:),'r-',tvec,sol(3,:),'g-',tvec,sol(4,:),'k-');
legend('W', 'Frq', 'Pcy', 'Pnu');
Periode = calculperiode(sol(1,:), tvec);
disp("Période calculée : ");
disp(Periode);

Sensibilite=sensitivity(par, tvec);
disp("Sensibilité calculée : ");
disp(Sensibilite);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function C = fonction_cout(p)
    
    // Conditions initiales
    x0 = [1; 0.1; 0.1; 0.1];
    t0 = 0;
    Nc = 30; // Nombre de cycles
    tvec = t0:0.05:Nc*24; // Temps de simulation de 0 à 30 cycles avec un pas de 0.05
    
    // Solution du système
    xsol = ode("stiff", x0, t0, tvec, list(clock_model_rhs, p));
    
    // Calcul de la période du système
    periodeS = calculperiode(xsol(1,:), tvec);
    
    // Calcul des amplitudes maximales et minimales
    Xmax = max(xsol(1,:));
    Xmin = min(xsol(1,:));
    
    // Comparaison avec les données observées
    observed_period = 24; // Période observée (à ajuster selon vos données)
    observed_Xmax = 14/16; // Amplitude maximale observée (à ajuster selon vos données)
    observed_Xmin = 1/16; // Amplitude minimale observée (à ajuster selon vos données)
    
    // Calcul de la fonction de coût
    C = (periodeS - observed_period)^2 + (Xmax - observed_Xmax)^2 + (Xmin - observed_Xmin)^2;
    
endfunction

// Données observées
observed_period = 24; // Période observée (à ajuster selon vos données)
observed_Xmax = 14/16; // Amplitude maximale observée (à ajuster selon vos données)
observed_Xmin = 1/16; // Amplitude minimale observée (à ajuster selon vos données)
data_obs = [observed_period; observed_Xmax; observed_Xmin];

// Paramètres initiaux
initial_guess_params = [a1_init, a2_init, a3_init, a4_init, K1_init, K4_init, gamma_init, gamma_m_init, gamma_P_init, n_init];

// Appel de la fonction d'optimisation pour minimiser J(p)
best_params = fminsearch(fonction_cout, initial_guess_params);

disp('Meilleurs paramètres :');
disp(best_params);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Question 6

data_FRQ=[2.5;5.5;9;5;4.3;0.5;3;7;5;3;1.2;0.9]; //à remplir à partir des graphiques présents sur le projet
data_frq=[9;4;2.5;0.5;0.4;9;3.8;1.5;1;1.2;0.8];
data_temps=[11;15;18;25;29;33;36;40;43;47;51];
xsol=ode("stiff",x0,t0,tvec,list(clock_model_rhs,par));
dt=0.01;
tvec=[0:dt:200];

//Données et modèles, amplitudes normalisées

Ptotal_mod=xsol(3,:)+xsol(4,:);
Pntotal=Ptotal_mod/max(Ptotal_mod);
data_FRQn=data_FRQ/max(data_FRQ);
data_frqn=data_frq/max(data_frq);
frq_mod=xsol(2,:)/max(xsol(2,:));

//Identifier "deux cycles" du modèle à partir d'un max de frq_mod
L=length(frq_mod);
i1=floor(0.25*L); //laisser passer le transcient
[mxfrq,imxfrq]=max(frq_mod(i1:L)); //imxfrq= indices dans la matrice xsol, chaque indice rpz la valeur d'une variable à un temps t
i2=imxfrq +48/dt;
tindexes=[imxfrq:1:i2];
figure(10);
offset=tvec(imxfrq)-data_temps(1);
plot(tvec(tindexes)+offset,frq_mod(tindexes),'K-');
plot(data_temps,data_frq,'K.');

//¨Partie principale 2

//Pour estimer les meilleurs paramètres

[petale, err]=fminsearch(fonction_cout(), par);

disp(petale);
disp(err);
//petale=vecteur optimisé, de même dimension que par()
//err=erreur finale de fct cout

// Création d'une fonction pour calculer la période en fonction de a3

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Création d'une fonction qui permet de calculer la période en fonction de a3
function [a3_values, periods] = compute_period_vs_a3(par, percentage_change)
    // Initialiser les valeurs
    a3 = par(3);
    initial_a3 = a3;
    a3_values = [initial_a3 * (1 - percentage_change), initial_a3, initial_a3 * (1 + percentage_change)];
    periods = zeros(1, length(a3_values));
    
    // Calculer la période pour chaque valeur de a3
    for i = 1:length(a3_values)
        par(3) = a3_values(i);
        sol = ode("stiff", x0, t0, tvec, list(clock_model_rhs, par));
        periods(i) = calculperiode(sol(1,:), tvec);
    end
endfunction

// Paramètres initiaux
initial_a3 = par(3);
percentage_change = 0.2; // on fait varier de ± 20% 

% Calculer la période pour différentes valeurs de a3
[a3_values, periods] = compute_period_vs_a3(par, percentage_change);

// Faire une figure de la période en fonction de a3 pour avoir une idée de l'impact
figure();
plot(a3_values, periods, 'bo-');
xlabel('a3');
ylabel('Période');
title('Variation de la période en fonction de a3');

