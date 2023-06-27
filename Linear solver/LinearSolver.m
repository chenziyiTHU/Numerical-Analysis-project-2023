% Linear Solver based on Lee and Greengard

% Each node structure contains: left&right endpoint, alpha_l, alpha_r, beta_l,
% beta_r, delta_l, delta_r, lambda_l, lambda_r, lambda
tic;

clear;
clc;
Initialization; % Initialization

TEST=1;
TOL=10^(-8);
C=4;
MaxLoop=6; % 最大循环上限
%DblFlag=false;

if q0==0
    s=(c-a)*zeta_l0*zeta_r0-zeta_l1*zeta_r0+zeta_r1*zeta_l0;
elseif q0==-1
    s=zeta_l1*zeta_r1*sinh(a-c)-zeta_l1*zeta_r0*cosh(a-c)+zeta_r1*zeta_l0*cosh(a-c)+zeta_l0*zeta_r0*sinh(c-a);
end

interval=struct('leftchild',0,'rightchild',0,'parent',0,'data',0,'exist',1);
interval.data=struct('Interval',[0,0],'alpha_l',0,'alpha_r',0, ...
    'beta_l',0,'beta_r',0,'delta_l',0,'delta_r',0, ...
    'lambda_l',0,'lambda_r',0,'lambda',0);

interval(1).data.Interval=[a,c]; 
interval(1).data.lambda_l=0;
interval(1).data.lambda_r=0;  
interval(1).data.lambda=1;

N=1;

R=0; % counting number of rounds

err=zeros(MaxLoop,1); % L2 error with the exact solution
err_infty=zeros(MaxLoop,1); % L_infty error with the exact solution

temp=0;
temp_discretization=0;

while(TEST>TOL)
    P_bar_inverse_f_tilde=zeros(K,N);
    P_bar_inverse_psi_l=zeros(K,N);
    P_bar_inverse_psi_r=zeros(K,N);

    % Step 1 (Local solver)
    M=0;
    for i=1:N
        if interval(i).leftchild==0 && interval(i).exist==1
            [alpha_l_Bi,alpha_r_Bi,...
                    beta_l_Bi,beta_r_Bi,...
                    delta_l_Bi,delta_r_Bi,...
                    P_bar_inverse_f_tilde_Bi,...
                    P_bar_inverse_psi_l_Bi,...
                    P_bar_inverse_psi_r_Bi]=Local_solver(interval(i).data.Interval(1),interval(i).data.Interval(2),K);

            interval(i).data.alpha_l=alpha_l_Bi;
            interval(i).data.alpha_r=alpha_r_Bi;
            interval(i).data.beta_l=beta_l_Bi;
            interval(i).data.beta_r=beta_r_Bi;
            interval(i).data.delta_l=delta_l_Bi;
            interval(i).data.delta_r=delta_r_Bi;

            P_bar_inverse_f_tilde(:,i)=P_bar_inverse_f_tilde_Bi;
            P_bar_inverse_psi_l(:,i)=P_bar_inverse_psi_l_Bi;
            P_bar_inverse_psi_r(:,i)=P_bar_inverse_psi_r_Bi;

            M=M+1; % number of leaf nodes
        end
    end

    % Step 2.A (Upward sweep)
    for i=N:-1:1
        if interval(i).leftchild>0 && interval(i).exist==1
            lc=interval(i).leftchild;
            rc=interval(i).rightchild;
            
            alpha_l_kl=interval(lc).data.alpha_l;
            alpha_l_kr=interval(rc).data.alpha_l;
            alpha_r_kl=interval(lc).data.alpha_r;
            alpha_r_kr=interval(rc).data.alpha_r;
            beta_l_kl=interval(lc).data.beta_l;
            beta_l_kr=interval(rc).data.beta_l;
            beta_r_kl=interval(lc).data.beta_r;
            beta_r_kr=interval(rc).data.beta_r;
            delta_l_kl=interval(lc).data.delta_l;
            delta_l_kr=interval(rc).data.delta_l;
            delta_r_kl=interval(lc).data.delta_r;
            delta_r_kr=interval(rc).data.delta_r;

            Delta=1-alpha_r_kr*beta_l_kl;
            alpha_l_k=(1-alpha_l_kr)*(alpha_l_kl-beta_l_kl*alpha_r_kr)/Delta+alpha_l_kr;
            alpha_r_k=alpha_r_kr*(1-beta_r_kl)*(1-alpha_l_kl)/Delta+alpha_r_kl;
            beta_l_k=beta_l_kl*(1-beta_r_kr)*(1-alpha_l_kr)/Delta+beta_l_kr;
            beta_r_k=(1-beta_r_kl)*(beta_r_kr-beta_l_kl*alpha_r_kr)/Delta+beta_r_kl;
            delta_l_k=(1-alpha_l_kr)/Delta*delta_l_kl+delta_l_kr+(alpha_l_kr-1)*alpha_l_kl/Delta*delta_r_kr;
            delta_r_k=(1-beta_r_kl)/Delta*delta_r_kr+delta_r_kl+(beta_r_kl-1)*alpha_r_kr/Delta*delta_l_kl;

            interval(i).data.alpha_l=alpha_l_k;
            interval(i).data.alpha_r=alpha_r_k;
            interval(i).data.beta_l=beta_l_k;
            interval(i).data.beta_r=beta_r_k;
            interval(i).data.delta_l=delta_l_k;
            interval(i).data.delta_r=delta_r_k;
        end
    end

    % Step 2.B (Downward sweep)
    for i=1:N
        lc=interval(i).leftchild;
        rc=interval(i).rightchild;
        if lc>0
            lambda_k=interval(i).data.lambda;
            lambda_l_k=interval(i).data.lambda_l;
            lambda_r_k=interval(i).data.lambda_r;

            lambda_kl=lambda_k;
            lambda_kr=lambda_k;
            lambda_l_kl=lambda_l_k;
            lambda_r_kr=lambda_r_k;

            alpha_r_kr=interval(rc).data.alpha_r;
            alpha_l_kl=interval(lc).data.alpha_l;
            beta_l_kl=interval(lc).data.beta_l;
            beta_r_kr=interval(rc).data.beta_r;
            delta_r_kr=interval(rc).data.delta_r;
            delta_l_kl=interval(lc).data.delta_l;

            A=[1, alpha_r_kr; beta_l_kl, 1];
            B=[lambda_r_k*(1-beta_r_kr)-lambda_k*delta_r_kr; lambda_l_k*(1-alpha_l_kl)-lambda_k*delta_l_kl];
            ss=A\B;
            lambda_r_kl=ss(1);
            lambda_l_kr=ss(2);

            interval(lc).data.lambda_l=lambda_l_kl;
            interval(lc).data.lambda_r=lambda_r_kl;
            interval(lc).data.lambda=lambda_kl;
            interval(rc).data.lambda_l=lambda_l_kr;
            interval(rc).data.lambda_r=lambda_r_kr;
            interval(rc).data.lambda=lambda_kr;
        end
    end

    % Step 2.C (Evaluation of solution of global integral equation)
    sigma=zeros(K,M);
    leaf_nodes=zeros(1,M);
    idx=0;
    for i=1:N
        if interval(i).leftchild==0 && interval(i).exist==1
            idx=idx+1;
            leaf_nodes(idx)=i;
        end
    end
    xx=zeros(M,1);
    for i=1:M
        xx(i)=interval(leaf_nodes(i)).data.Interval(1);
    end
    [~,index]=sort(xx);
    for idx=1:M
        i=leaf_nodes(index(idx));
        l=interval(i).data.Interval(1);
        r=interval(i).data.Interval(2);
        tau=zeros(K,1);
        for j=1:K
            tau(j)=(r-l)/2*cos((2*K-2*j+1)*pi/(2*K))+(l+r)/2;
        end
        for j=1:K
            sigma(j,idx)=P_bar_inverse_f_tilde(j,i)+interval(i).data.lambda_l*P_bar_inverse_psi_l(j,i)+interval(i).data.lambda_r*P_bar_inverse_psi_r(j,i);
        end
    end
    % 注意这里的区间B_i对应的是结构体interval(leaf_nodes(index(i)))

    % Step 3.A (Evaluation of solution)
    J_l=zeros(M,1);
    J_r=zeros(M,1);
    J_l(1)=0;
    J_r(M)=0;
    for i=1:M-1
        J_l(i+1)=J_l(i)+interval(leaf_nodes(index(i))).data.delta_l+interval(leaf_nodes(index(i))).data.lambda_l*interval(leaf_nodes(index(i))).data.alpha_l+interval(leaf_nodes(index(i))).data.lambda_r*interval(leaf_nodes(index(i))).data.beta_l;
    end
    for i=M:-1:2
        J_r(i-1)=J_r(i)+interval(leaf_nodes(index(i))).data.delta_r+interval(leaf_nodes(index(i))).data.lambda_l*interval(leaf_nodes(index(i))).data.alpha_r+interval(leaf_nodes(index(i))).data.lambda_r*interval(leaf_nodes(index(i))).data.beta_r;
    end

    % Step 3.B (Evaluation of solution)
    u=zeros(K,M);
    discretization=zeros(K,M);
    for i=1:M
        l=interval(leaf_nodes(index(i))).data.Interval(1);
        r=interval(leaf_nodes(index(i))).data.Interval(2);
        tau=zeros(K,1);
        for j=1:K
            tau(j)=(r-l)/2*cos((2*K-2*j+1)*pi/(2*K))+(l+r)/2;
        end
        aa=g_l(tau(1:K)).*sigma(1:K,i);
        bb=g_r(tau(1:K)).*sigma(1:K,i);
        for j=1:K
            u(j,i)=u_i(tau(j))+g_r(tau(j))/s*(J_l(i)+LIntegral(l,r,tau(j),aa))+g_l(tau(j))/s*(RIntegral(l,r,tau(j),bb)+J_r(i));
            discretization(j,i)=tau(j);
        end
    end

    % Step 4 (Mesh refinement)
    if R==0
        TEST=1;
    else
        TEST=Test(u,discretization,temp,temp_discretization);
    end
    temp=u;
    temp_discretization=discretization;
    S=zeros(M,1);
    for i=1:M
        chebcoef_sigma1=0;
        chebcoef_sigma2=0;
        chebcoef_sigma3=0;
        for j=1:K
            chebcoef_sigma1=chebcoef_sigma1+2/K*sigma(j,i)*cos((K-1)*acos(discretization(j,i)));
            chebcoef_sigma2=chebcoef_sigma2+2/K*sigma(j,i)*cos((K-2)*acos(discretization(j,i)));
            chebcoef_sigma3=chebcoef_sigma3+2/K*sigma(j,i)*cos((K-3)*acos(discretization(j,i)));
        end
        S(i)=abs(chebcoef_sigma2)+abs(chebcoef_sigma1-chebcoef_sigma3);
    end
    S_div=max(S)/(2^C);
    for i=1:M
        if S(i)>=S_div
            l=interval(leaf_nodes(index(i))).data.Interval(1);
            r=interval(leaf_nodes(index(i))).data.Interval(2);
            interval(N+1)=struct('leftchild',0,'rightchild',0,'parent',leaf_nodes(index(i)),'data',0,'exist',1);
            interval(N+1).data=struct('Interval',[l,(l+r)/2],'alpha_l',0,'alpha_r',0, ...
                    'beta_l',0,'beta_r',0,'delta_l',0,'delta_r',0, ...
                    'lambda_l',0,'lambda_r',0,'lambda',0);
            interval(leaf_nodes(index(i))).leftchild=N+1;
            interval(N+2)=struct('leftchild',0,'rightchild',0,'parent',leaf_nodes(index(i)),'data',0,'exist',1);
            interval(N+2).data=struct('Interval',[(l+r)/2,r],'alpha_l',0,'alpha_r',0, ...
                    'beta_l',0,'beta_r',0,'delta_l',0,'delta_r',0, ...
                    'lambda_l',0,'lambda_r',0,'lambda',0);
            interval(leaf_nodes(index(i))).rightchild=N+2;
            N=N+2;
        elseif i<M && interval(leaf_nodes(index(i))).parent==interval(leaf_nodes(index(i+1))).parent && S(i)+S(i+1)<S_div/(2^K) && interval(leaf_nodes(index(i))).parent>0
            interval(leaf_nodes(index(i))).exist=0;
            interval(leaf_nodes(index(i+1))).exist=0;
            m=interval(leaf_nodes(index(i))).parent;
            interval(m).leftchild=0;
            interval(m).rightchild=0;
        end
    end

    R=R+1;

    % Compute the L2 error with the exact solution
    I=0;
    for i=1:M
        l=interval(leaf_nodes(index(i))).data.Interval(1);
        r=interval(leaf_nodes(index(i))).data.Interval(2);
        err(R)=err(R)+Integral(l,r,(exactsolution(discretization(:,i))-u(:,i)).^2);
    end
    for i=1:M
        I=I+Integral(l,r,(exactsolution(discretization(:,i))).^2);
    end
    err(R)=sqrt(err(R)/I);

    % L_infty error
    for i=1:M
        for j=1:K
            err_infty(R)=max(err_infty(R),abs(exactsolution(discretization(j,i))-u(j,i)));
        end
    end
    
    if R==MaxLoop
        break
    end

end

toc

figure;
plot(discretization(:),u(:))
ylim([-2,2]);
figure;
semilogy(1:R,err(1:R),'O-')
ylim([10^(-10),10^10]);
xlabel('step','FontSize',15)
ylabel('L^2 Error','FontSize',15)
figure;
semilogy(1:R,err_infty(1:R),'O-')
ylim([10^(-10),10^10]);
xlabel('step','FontSize',15)
ylabel('L^{\infty} Error','FontSize',15)

