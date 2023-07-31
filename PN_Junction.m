GIUNZIONE P - N
Na>Nd, Condizioni di equilibrio termodinamico.
clc; clear all; close all;
%Dati e costanti
Nc=2.8e19; %cm^-3
Nv=1.04e19; %cm^-3
NA=2e16; %cm^-3
ND=1e16; %cm^-3
q=1.602e-19; %C
epsilon0=8.85e-12; %F/m
epsilonS=11.9*epsilon0;
T=300;%K
K=8.6167e-5; %eV/K
K1=K*q; %J/K
ni=1.45e10; %cm^-3
Vt=K1*T/q; %V - tensione termica

%Relazione di Boltzmann 
Vbi=K1*T*(1/q)*log(NA*ND/ni^2); %V
VBI_EV_C=Vbi*1/q; %eV/C
xn=10^-3*sqrt(  2 * epsilonS * (1/q) * Vbi * 1/(ND) * ( 1 / ( ( ND/NA ) + 1 ) )  );
xp=10^-3*sqrt(  2 * epsilonS * (1/q) * Vbi * 1/(NA) * ( 1 / ( ( NA/ND ) + 1 ) )  );
W=10^-3*sqrt(2*epsilonS*(1/q)*(Vbi)*(1/NA + 1/ND));
fprintf("Potenziale di built-in Vbi = %f V ",Vbi)
fprintf("Spessore della SCR nella regione di tipo N = %f mm = %f µm",xn*1000,xn*1000*1000)
fprintf("Spessore della SCR nella regione di tipo P = %f mm = %f µm",xp*1000,xp*1000*1000)
fprintf("Spessore della SCR = %f mm = %f µm",W*1000,W*1000*1000)

%Densità di carica volumetrica 
x1=linspace(-xp,0,50);
x2=linspace(0,xn,50);
rho_p=-q*NA+zeros(1,length(x1)); %C/cm^3
rho_n=q*ND+zeros(1,length(x2));  %C/cm^3
figure(1)
area(x1,rho_p);
hold on
grid on 
area(x2,rho_n);
xline(0,'k',"LineWidth",1.2)
yline(0,'w')
axis padded
ylabel("Densità di carica volumetrica ρ(x) [C/cm^3]",'FontWeight','bold')
xlabel("Asse x [m] ",'FontWeight','bold')
legend("q*N_A","q*N_D","x_j",'Location','best')

%Campo elettrico V/cm
E_p=-100*100*q*NA*(x1+xp)/epsilonS+zeros(1,length(x1));
E_n=100*100*q*ND*(x2-xn)/epsilonS+zeros(1,length(x2));
fprintf("|Emax|= %f V/cm\n",q*ND*xn*100*100/epsilonS)
figure(2)
plot(x1,E_p,"LineWidth",1)
hold on
plot(x2,E_n,"LineWidth",1)
xline(0,'k',"LineWidth",1.2)
axis padded
grid on 
ylabel("Campo elettrico ℰ(x) [V/cm]",'FontWeight','bold')
xlabel("Asse x [m] ",'FontWeight','bold')
legend("-qN_A(x+xp)/ε_s","qN_D(x-xn)/ε_s","Location","best")
%L'area sottesa dal campo elettrico cambiata di segno deve essere uguale
%alla Vbi
Area=100*W*max(abs(E_p))/2 %V
hold off

%Potenziale interno Vbi=Psi(xn)=Psi_p+Psi_n
Psi_p=100*100*100*q*NA*1/(2*epsilonS)*(x1+xp).^2; %V
Psi_n=100*100*100*( (-q)*ND*1/(2*epsilonS)*x2.^2 +q*ND*1/(epsilonS)*xn.*x2+q*NA*1/(2*epsilonS)*xp^2); %V
figure(3)
plot(x1,Psi_p,"LineWidth",1)
hold on
plot(x2,Psi_n,"LineWidth",1)
xline(0,'k',"LineWidth",1.2)
grid on  
axis padded
ylabel("Potenziale interno Ψ(x) [V]",'FontWeight','bold')
xlabel("Asse x [m] ",'FontWeight','bold')
legend("Ψ_p(x)","Ψ_n(x)","x_j","Location","best")
hold off

%Diagramma a bande di energia in eV
EG=1.12; %eV
Ec0=EG/2;
Ev0=-EG/2;
PSI_x=[Psi_p Psi_n];
x=[x1 x2];
Ec=Ec0-PSI_x;
Ev=Ev0-PSI_x;
Ei0=(Ec0+Ev0)/2 - K*T/2 * log(Nc/Nv);
Ei=Ei0-PSI_x;
figure(4)
plot(x,Ec,x,Ev,'Color','b')
hold on
plot(x,Ei,'--m')
grid on 
axis padded
ylabel("[eV]",'FontWeight','bold')
xlabel("Asse x [m] ",'FontWeight','bold')
legend("E_c","E_v","E_i",'Location','best')

%Concentrazione dei portatori di carica nella SCR
p=NA*exp(-PSI_x/Vt);
n=ni^2/NA * exp(PSI_x/Vt);
figure(5)
plot(x,p,'r',x,n,'b')
grid on 
hold on 
axis padded
x_np=linspace(-2*xp,-xp,50);
x_pn=linspace(xn,xn+xp,50);
%A conferma della legge della giunzione vediamo che non ci sono nè accumuli
%nè asportazioni di minoritari ai bordi della SCR
pn_xn=p(end) %pn(xn)=pno
pno=ni^2/ND
np_xp=n(1)   %np(-xp)=npo
npo=ni^2/NA
NA_v=NA+zeros(1,length(x_np));
ND_v=ND+zeros(1,length(x_pn));
plot(x_np,NA_v,'r:',x_pn,ND_v,'b:');
plot(x_pn,ni^2./NA_v,'r--',x_np,ni^2./ND_v,'b--')
xline(xn,'k')
xline(0,'k','LineWidth',1.3)
xline(-xp,'k')
xlabel("Asse x [m] ",'FontWeight','bold')
ylabel("Concentrazione di elettroni liberi e lacune [cm^-3]",'FontWeight','bold')
legend("p(x)","n(x)","p_{po} = N_A","n_{no} = N_D","p_{no} = n_i^2/N_D ","n_{po} = n_i^2/N_A","(-xp) ed xn ","x_j",'Location','best')
hold off 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6)
semilogy(x,p,'r',x,n,'b')
grid on 
hold on
axis padded
title("Andamento di n e p all'interno della SCR in scala logaritmica","FontWeight","bold")
xlabel("Asse x [m] ",'FontWeight','bold')
ylabel("Concentrazioni [cm^{-3}]",'FontWeight','bold')
yline(npo,'b--',"LineWidth",1)
yline(pno,'r:',"LineWidth",1.3)
yline(NA,'r:',"LineWidth",1.3)
yline(ND,'b--',"LineWidth",1)
legend("p(x)","n(x)","n_{po} = n_i^2/N_A","p_{no} = n_i^2/N_D ","p_{po} = N_A","n_{no} = N_D",'Location','east')
hold off

%Dipendenze dai principali parametri 
%Vbi=f(T,NA);
T_v=300:50:500;
ni_T=ni_function(T_v);
ND_v2=[ND/10 ND ND*10 ND*100];
for i=1:length(ND_v2)
Vbi_T(i,:)=K1.*T_v.*(1/q).*log(NA.*ND_v2(i)./(ni_T.^2));
end
figure(7)
plot(T_v,Vbi_T)
grid on 
hold on 
xlabel("Temperatura [K] ","FontWeight","bold")
ylabel("Vbi [V] ","FontWeight","bold")
legend("N_D=1e15","N_D=1e16","N_D=1e17","N_D=1e18",'Location', "best")
hold off

%W=f(NA,T)
NA_v2=logspace(15,19,30);
Vbi_NA_T=VBI_function(NA_v2,T_v);
ND_v1=1e15;
for i=1:length(NA_v2)
W_NA(i,:)=sqrt(2*epsilonS*(1/q).*Vbi_NA_T(i,:).*(1./NA_v2(i) + 1/ND))/1000; %metri
W_NA_ND_v1(i,:)=sqrt(2*epsilonS*(1/q).*Vbi_NA_T(i,:).*(1./NA_v2(i) + 1/ND_v1))/1000; %metri
end
figure(8)
for i=1:length(T_v)
semilogx(NA_v2,W_NA(:,i))
hold on
grid on 
end
xlabel(" N_A [cm^{-3}] ","FontWeight","bold")
ylabel(" W [m] ","FontWeight","bold")
legend("T_1=300 K","T_2=350 K","T_3=400 K","T_4=450 K","T_5=500 K",'Location', "best")
