function error=Test(u,u_discretization,v,v_discretization)
x1=u_discretization(:);
x2=v_discretization(:);
u=u(:);
v=v(:);
l1=length(x1);
l2=length(x2);
diff=0;
sum=0;
for i=1:l2
    for j=1:l1-1
        if x2(i)>=x1(j) && x2(i)<=x1(j+1)
            t=(u(j+1)*(x2(i)-x1(j))+u(j)*(x1(j+1)-x2(i)))/(x1(j+1)-x1(j));
            diff=max(diff,abs(t-v(i)));
            sum=max(sum,abs(t+v(i)));
            break
        end
    end
end
error=diff/sum;