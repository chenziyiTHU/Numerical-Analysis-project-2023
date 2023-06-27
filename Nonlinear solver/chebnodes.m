function x=chebnodes(l,r,K)
x=zeros(K,1);
for i=1:K
    x(i)=(r-l)/2*cos((2*K-2*i+1)*pi/(2*K))+(r+l)/2;
end