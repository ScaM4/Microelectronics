clc;
clear;
close all;
% Funzione di distribuzione di probabilità di Fermi - Dirac
% F(E)=1/(1+exp(E-Ef/KT)) -> Probabilità che un elettrone occupi uno stato
% quanto al livello energetico E

K=8.617e-5;             %eV/K
T=[0.001,100:200:1500]; %K 
x=-0.5:0.01:0.5;        %E-Ef eV

F=zeros(length(T),length(x));
for i=1:length(T)
 F(i,:)=1./(exp(x/(K*T(i)))+1);
 plot(x,F(i,:),'LineWidth',1)
 grid on
 hold on
 title([],'FontWeight','bold');
 xlabel('E-Ef [eV]','FontWeight','bold')
 ylabel('F(E)','FontWeight','bold')
end
legend('T=0.001K','T=100K','T=300K','T=500K','T=700K','T=900K','T=1100K','T=1300K','T=1500K')
annotation('arrow',[0.5059 0.3484],[0.915 0.7375])
annotation('arrow',[0.5236 0.689],[0.1119 0.2966])
hold off
