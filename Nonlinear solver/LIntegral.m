function L=LIntegral(l,r,x,values)
K=length(values);
nodes=zeros(K,1);
alpha=zeros(K,1);
for j=1:K
    nodes(j)=cos((2*K-2*j+1)*pi/(2*K));
end
for k=2:K
    for j=1:K
        alpha(k)=alpha(k)+2/K*values(j)*cos((k-1)*(2*K-2*j+1)*pi/(2*K));
    end
end
for j=1:K
        alpha(1)=alpha(1)+1/K*values(j);
end

a=zeros(K+1,1);
for k=3:K-1
    a(k)=1/(2*(k-1))*(alpha(k-1)-alpha(k+1));
end
a(K)=1/(2*(K-1))*alpha(K-1);
a(K+1)=1/(2*K)*alpha(K);
a(2)=1/2*(2*alpha(1)-alpha(3));
for k=2:K+1
    a(1)=a(1)+(-1)^k*a(k);
end

x0=2*x/(r-l)-(r+l)/(r-l);
L=0;
for k=1:K+1
    L=L+a(k)*cos((k-1)*acos(x0));
end
L=L*(r-l)/2;