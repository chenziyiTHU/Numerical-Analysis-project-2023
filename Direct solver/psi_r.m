function val=psi_r(x)
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
    s=(c-a)*zeta_l0*zeta_r0-zeta_l1*zeta_r0+zeta_r1*zeta_l0;
    val=(p_tilde(x)*zeta_l0+ q_tilde(x)*g_l(x))/s;
elseif q0==-1
    s=(zeta_l1*zeta_r1-zeta_l0*zeta_r0)*sinh(a-c)+(zeta_l0*zeta_r1-zeta_l1*zeta_r0)*cosh(a-c);
    val=(p_tilde(x)*(zeta_l1*sinh(x-c)-zeta_l0*cosh(x-c))-q_tilde(x)*(zeta_l1*cosh(x-c)-zeta_l0*sinh(x-c)))/s;
end
