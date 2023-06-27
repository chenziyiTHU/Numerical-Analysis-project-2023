function R=RIntegral(l,r,x,values)
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

b=zeros(K+1,1);
for k=3:K-1
    b(k)=1/(2*(k-1))*(alpha(k+1)-alpha(k-1));
end
b(K)=-1/(2*(K-1))*alpha(K-1);
b(K+1)=-1/(2*K)*alpha(K);
b(2)=1/2*(alpha(3)-2*alpha(1));
for k=2:K+1
    b(1)=b(1)-b(k);
end

x0=2*x/(r-l)-(r+l)/(r-l);
R=0;
for k=1:K+1
    R=R+b(k)*cos((k-1)*acos(x0));
end
R=R*(r-l)/2;