% Direct solver
% Still consider euquation of the form: u''(x)+p(x)u'(x)+q(x)u(x)=f
clear;
clc;

Initialization;
if q0==0
    s=(c-a)*zeta_l0*zeta_r0-zeta_l1*zeta_r0+zeta_r1*zeta_l0;
elseif q0==-1
    s=zeta_l1*zeta_r1*sinh(a-c)-zeta_l1*zeta_r0*cosh(a-c)+zeta_r1*zeta_l0*cosh(a-c)+zeta_l0*zeta_r0*sinh(c-a);
end

% Partition of the interval
n=35; % Number of partitions
aa=zeros(n,1);
cc=zeros(n,1); % Require cc(i)=aa(i+1)

% For convenience, consider equidistant intervals
for i=1:n
    aa(i)=(c-a)/n*(i-1)+a;
    cc(i)=(c-a)/n*i+a;
end

x=zeros(K*n,1);
PSI_l=zeros(K*n);
PSI_r=zeros(K*n);
for i=1:n
    tt=chebnodes(aa(i),cc(i),K);
    for j=1:K
        x((i-1)*K+j)=tt(j);
        PSI_l((i-1)*K+j,(i-1)*K+j)=psi_l(tt(j));
        PSI_r((i-1)*K+j,(i-1)*K+j)=psi_r(tt(j));
    end
end % Discretization

e=eye(K*n);

% Calculate the matrix P
P=zeros(K*n);

for i=1:n
    for j=1:K
        % Consider P act on vector e_{(i-1)*K+j}
        y=e(:,(i-1)*K+j);
        value=zeros(K*n,1);
        for k=1:n
            for l=1:K
                value((k-1)*K+l)=y((k-1)*K+l);
                if k<i
                    value((k-1)*K+l)=value((k-1)*K+l)...
                        +PSI_r((k-1)*K+l,(k-1)*K+l)*Integral(aa(i),cc(i),g_r(x(((i-1)*K+1):(i*K))).*y(((i-1)*K+1):(i*K)));
                elseif k==i
                    value((k-1)*K+l)=value((k-1)*K+l)...
                        +PSI_l((k-1)*K+l,(k-1)*K+l) * LIntegral(aa(k), cc(k), x((k-1)*K+l), g_l(x(((k-1)*K+1):(k*K))).*y(((k-1)*K+1):(k*K)))...
                        +PSI_r((k-1)*K+l,(k-1)*K+l) * RIntegral(aa(k), cc(k), x((k-1)*K+l), g_r(x(((k-1)*K+1):(k*K))).*y(((k-1)*K+1):(k*K)));
                else
                    value((k-1)*K+l)=value((k-1)*K+l)...
                        +PSI_l((k-1)*K+l,(k-1)*K+l)*Integral(aa(i),cc(i),g_l(x(((i-1)*K+1):(i*K))).*y(((i-1)*K+1):(i*K)));
                end
            end
        end

        P(:,(i-1)*K+j)=value;
    end
end

% Compute sigma
b=f_tilde(x);
sigma=P\b;

% Compute J_l and J_r
J_l=zeros(n);
J_r=zeros(n);
J_l(1)=0;
J_r(n)=0;
for i=1:n-1
    J_l(i+1)=J_l(i)+Integral(aa(i),cc(i),g_l(x(((i-1)*K+1):(i*K))).*sigma(((i-1)*K+1):(i*K)));
end
for i=n:-1:2
    J_r(i-1)=J_r(i)+Integral(aa(i),cc(i),g_r(x(((i-1)*K+1):(i*K))).*sigma(((i-1)*K+1):(i*K)));
end

% Compute u
u=u_i(x);
for i=1:n
    for j=1:K
        u((i-1)*K+j)=u((i-1)*K+j)+g_r(x((i-1)*K+j))/s*(J_l(i)+LIntegral(aa(i),cc(i),x((i-1)*K+j),g_l(x(((i-1)*K+1):(i*K))).*sigma(((i-1)*K+1):(i*K))))...
            +g_l(x((i-1)*K+j))/s*(J_r(i)+RIntegral(aa(i),cc(i),x((i-1)*K+j),g_r(x(((i-1)*K+1):(i*K))).*sigma(((i-1)*K+1):(i*K))));
    end
end

plot(x,u)




