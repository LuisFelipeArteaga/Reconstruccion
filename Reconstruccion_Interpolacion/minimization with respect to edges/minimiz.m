function [f]= minimiz(k,n)

f=zeros(1,2*n-1);
for i=1:1:(n-1)
    f(2*i-1)=(k(i)^2+k(i+1)^2)/2;
    f(2*i)=k(i+1)^2;
end

f(2*n-1)=(k(n)^2+k(1)^2)/2;
