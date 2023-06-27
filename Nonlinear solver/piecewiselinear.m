function val=piecewiselinear(x,u,x0)
K=lengt(x);
for j=1:K-1
    if x0>=x(j) && x0<x(j+1)
        val=(u(j+1)*(x0-x(j))+u(j)*(x(j+1)-x0))/(x(j+1)-x(j));
        break
    end
end