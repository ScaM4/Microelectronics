function ni_T=ni_function(T)
ni=1.45e10; %cm^-3
EG=1.206; %eV
K=8.6167e-5; %eV/K
ni_T=ni.*(T/300).^(1.5).*exp(-(EG/(2*K))*(1./T - 1/300));
end