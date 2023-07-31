function Vbi_NA_T=VBI_function(NA_vettore,T_v)
K=8.6167e-5; %eV/K
q=1.602e-19; %C
K1=K*q; %J/K
ND=1e16; %cm^-3
for i=1:length(NA_vettore)
Vbi_NA_T(i,:)=K1.*T_v.*(1/q).*log(NA_vettore(i).*ND./(ni_function(T_v).^2));
end 
end