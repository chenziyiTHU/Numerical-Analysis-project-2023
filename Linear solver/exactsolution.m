function u=exactsolution(x)
epsilon=evalin('base', 'epsilon');
u=1/erf(1/sqrt(epsilon))*erf(x/sqrt(epsilon));