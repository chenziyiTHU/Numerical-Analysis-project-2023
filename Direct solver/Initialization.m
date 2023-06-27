% Initialization
% Constants of the equation:
% zeta_l0, zeta_l1, zeta_r0, zeta_r1, Gamma_l, Gamma_r;
% Functions: p,q,f

% Computed functions: sigma(x), psi_l(x), psi_r(x), g_l(x), g_r(x).

a=-1;
c=1;
zeta_l0=1;
zeta_l1=0;
zeta_r0=1;
zeta_r1=0;
Gamma_l=-1;
Gamma_r=1;
q0=0;
K=16; % Number of chebyshev nodes on each interval
epsilon=10^(-5);