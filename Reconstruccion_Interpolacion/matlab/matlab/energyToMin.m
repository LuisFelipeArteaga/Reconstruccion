function e=energyToMin(arg,K,C)

x=arg(1:(end-1)).*C
[n,~]=size(x);

p=2;

e=1/2*norm(x-K,p)^p-arg(end)*(ones(1,n)*x-(n-2)*180);

return