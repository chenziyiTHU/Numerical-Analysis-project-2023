function val=f_tilde(x)
a=evalin('base', 'a');
c=evalin('base', 'c');
zeta_l0=evalin('base', 'zeta_l0');
zeta_r0=evalin('base', 'zeta_r0');
zeta_l1=evalin('base', 'zeta_l1');
zeta_r1=evalin('base', 'zeta_r1');
Gamma_l=evalin('base', 'Gamma_l');
Gamma_r=evalin('base', 'Gamma_r');
m=(zeta_r0*Gamma_l-zeta_l0*Gamma_r)/((a*zeta_l0+zeta_l1)*zeta_r0-(c*zeta_r0+zeta_r1)*zeta_l0);
n=((a*zeta_l0+zeta_l1)*Gamma_r-(c*zeta_r0+zeta_r1)*Gamma_l)/((a*zeta_l0+zeta_l1)*zeta_r0-(c*zeta_r0+zeta_r1)*zeta_l0);
val=f(x)-(p(x)*m+q(x)*(m*x+n));