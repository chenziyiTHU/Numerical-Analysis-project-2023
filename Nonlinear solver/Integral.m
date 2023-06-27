function I=Integral(l,r,values)
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
I=0;
for k=1:K
    if rem(k,2)~=0
        I=I+2*alpha(k)/(1-(k-1)^2);
    end
end
I=I*(r-l)/2;