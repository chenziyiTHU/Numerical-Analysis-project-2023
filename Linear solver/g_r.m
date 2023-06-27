function val=g_r(x)
a=evalin('base', 'a');
c=evalin('base', 'c');
zeta_l0=evalin('base', 'zeta_l0');
zeta_r0=evalin('base', 'zeta_r0');
zeta_l1=evalin('base', 'zeta_l1');
zeta_r1=evalin('base', 'zeta_r1');
Gamma_l=evalin('base', 'Gamma_l');
Gamma_r=evalin('base', 'Gamma_r');
q0=evalin('base', 'q0');
if q0==0
    val=zeta_r0*(x-c)-zeta_r1;
elseif q0==-1
    val=zeta_r1*cosh(x-c)-zeta_r0*sinh(x-c);
end