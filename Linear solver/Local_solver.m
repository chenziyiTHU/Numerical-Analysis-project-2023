% Step 1: Local solver
% l: left endpoint;
% r: right endpoint;
% K: order of accuracy.
function [alpha_l,alpha_r,...
    beta_l,beta_r,...
    delta_l,delta_r,...
    P_bar_inverse_f_tilde,...
    P_bar_inverse_psi_l,...
    P_bar_inverse_psi_r]=Local_solver(l,r,K)


% (1) Determine the locations of the scaled Chebyshev nodes.
tau=zeros(K,1);
for j=1:K
    tau(j)=(r-l)/2*cos((2*K-2*j+1)*pi/(2*K))+(l+r)/2;
end

% (2) Evaluate p_tilde, q_tilde, f_tilde at the scaled Chebyshev nodes.
p_tilde_chebnodes=zeros(K,1);
q_tilde_chebnodes=zeros(K,1);
f_tilde_chebnodes=zeros(K,1);
for j=1:K
    p_tilde_chebnodes(j)=p_tilde(tau(j));
    q_tilde_chebnodes(j)=q_tilde(tau(j));
    f_tilde_chebnodes(j)=f_tilde(tau(j));
end

% (3) Construct the linear system P_bar
psi_l_chebnodes=zeros(K,1);
psi_r_chebnodes=zeros(K,1);
for j=1:K
    psi_l_chebnodes(j)=psi_l(tau(j));
    psi_r_chebnodes(j)=psi_r(tau(j));
end
g_l_chebnodes=zeros(K,1);
g_r_chebnodes=zeros(K,1);
for j=1:K
    g_l_chebnodes(j)=g_l(tau(j));
    g_r_chebnodes(j)=g_r(tau(j));
end
P_bar=zeros(K);
for j=1:K
    e=zeros(K,1);
    e(j)=1;
    %P_bar(:,j)=e+psi_l_chebnodes(j)*LIntegral(l,r,tau(1:K),g_l_chebnodes.*e)+psi_r_chebnodes(j)*RIntegral(l,r,tau(1:K),g_r_chebnodes.*e);
    value=zeros(K,1);
    for ll=1:K
        value(ll)=e(ll)+psi_l_chebnodes(ll)*LIntegral(l,r,tau(ll),g_l_chebnodes.*e)+psi_r_chebnodes(ll)*RIntegral(l,r,tau(ll),g_r_chebnodes.*e);
    end
    P_bar(:,j)=value;
end

% (4) Evaluate P_bar^(-1)f_tilde, P_bar^(-1)psi_l, P_bar^(-1)psi_r.
P_bar_inverse_f_tilde=P_bar\f_tilde_chebnodes;
P_bar_inverse_psi_l=P_bar\psi_l_chebnodes;
P_bar_inverse_psi_r=P_bar\psi_r_chebnodes;

% (5) Evaluate alpha_l, alpha_r, beta_l, beta_r, delta_l, delta_r.
alpha_l=Integral(l,r,g_l_chebnodes.*P_bar_inverse_psi_l);
alpha_r=Integral(l,r,g_r_chebnodes.*P_bar_inverse_psi_l);
beta_l=Integral(l,r,g_l_chebnodes.*P_bar_inverse_psi_r);
beta_r=Integral(l,r,g_r_chebnodes.*P_bar_inverse_psi_r);
delta_l=Integral(l,r,g_l_chebnodes.*P_bar_inverse_f_tilde);
delta_r=Integral(l,r,g_r_chebnodes.*P_bar_inverse_f_tilde);




