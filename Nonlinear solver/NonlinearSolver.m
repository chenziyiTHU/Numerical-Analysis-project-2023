% Nonlinear solver
% Consider the example of lambert W function
tic;

clear;
clc;

a=exp(1);
c=2*exp(2);
times=7; % times of iteration

u0=struct('u',@(x) 0,'du',@(x) 0,'ddu',@(x) 0);

u0(1).u=@(x) (x+2*exp(2)-2*exp(1))/(2*exp(2)-exp(1));
u0(1).du=@(x) 1/(2*exp(2)-exp(1));
u0(1).ddu=@(x) 0; % Step 1: initial guess

f=@(x,y,z) z/(x*(y+1)^2)-y/(x^2*(y+1));
f2=@(x,y,z) -y/(x^2*(y+1)^2)-2*z/(x*(y+1)^3);
f3=@(x,y,z) 1/(x*(y+1)^2);
s=2*exp(2)-exp(1);

err1=zeros(1,times);
X=exp(1):(2*exp(2)-exp(1))/100:2*exp(2);
y=@(x) lambertw(x);

for RR=1:times
    p_tilde=@(x) -f3(x,u0(RR).u(x),u0(RR).du(x));
    p=p_tilde;
    q_tilde=@(x) -f2(x,u0(RR).u(x),u0(RR).du(x));
    q=q_tilde;
    f_tilde=@(x) -u0(RR).ddu(x)+f(x,u0(RR).u(x),u0(RR).du(x));
    g_l=@(x) x-exp(1);
    g_r=@(x) x-2*exp(2);
    psi_l=@(x) (p_tilde(x)+q_tilde(x)*g_r(x))/s;
    psi_r=@(x) (p_tilde(x)+q_tilde(x)*g_l(x))/s;
    u_i=@(x) 0;

    LinearSolver; % Step 2: solving linear ordinary differential equation
    
    le=zeros(M,1);
    re=zeros(M,1);
    for i=1:M
        le(i)=interval(leaf_nodes(index(i))).data.Interval(1);
        re(i)=interval(leaf_nodes(index(i))).data.Interval(2);
    end
    
    % Step 3: iteration
    ddu_solution=@(x) f3(x,u0(RR).u(x),u0(RR).du(x))*du_solution(x,J_l,J_r,le,re,sigma,K,g_l,g_r,s)+f2(x,u0(RR).u(x),u0(RR).du(x))*u_solution(x,J_l,J_r,le,re,sigma,K,g_l,g_r,s)+f_tilde(x);
    u0(RR+1).u=@(x) u0(RR).u(x)+u_solution(x,J_l,J_r,le,re,sigma,K,g_l,g_r,s);
    u0(RR+1).du=@(x) u0(RR).du(x)+du_solution(x,J_l,J_r,le,re,sigma,K,g_l,g_r,s);
    u0(RR+1).ddu=@(x) u0(RR).ddu(x)+ddu_solution(x);
    
    for i=1:101
        err1(RR)=max(err1(RR),abs(u0(RR).u(X(i))-y(X(i))));
    end
end

Y=zeros(101,1);
for i=1:101
    Y(i)=u0(times+1).u(X(i));
end

figure;
plot(X,y(X),'-r')
hold on
plot(X,Y,'-b')

figure;
semilogy(1:times,err1,'o-b');
xlabel('iteration step','FontSize',15)
ylabel('L^\infty Error','FontSize',15)
% hold on
% semilogy(1:times,exp(-4*(1:times)),'o-r');

toc

function val=u_solution(x,J_l,J_r,le,re,sigma,K,g_l,g_r,s)
    val=0;
    M=length(le);
    for i=1:M
        if x>le(i) && x<=re(i)
            val=g_r(x)/s.*(J_l(i)+LIntegral(le(i),re(i),x,g_l(chebnodes(le(i),re(i),K)).*sigma(1:K,i)))...
                +g_l(x)/s.*(RIntegral(le(i),re(i),x,g_r(chebnodes(le(i),re(i),K)).*sigma(1:K,i))+J_r(i));
            break
        end
    end
end

function val=du_solution(x,J_l,J_r,le,re,sigma,K,g_l,g_r,s)
    val=0;
    M=length(le);
    for i=1:M
        if x>le(i) && x<=re(i)
            val=1/s.*(J_l(i)+LIntegral(le(i),re(i),x,g_l(chebnodes(le(i),re(i),K)).*sigma(1:K,i)))...
                +1/s.*(RIntegral(le(i),re(i),x,g_r(chebnodes(le(i),re(i),K)).*sigma(1:K,i))+J_r(i));
            break
        end
    end
end





