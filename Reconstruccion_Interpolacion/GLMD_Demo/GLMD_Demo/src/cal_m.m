%======================================================
% Calculate Corresponding Matrix
%======================================================
function [M] = cal_m(source,Nm,target,Nf,T) 
[n, dims] = size (source);
[m, dims] = size (target);
M0 = ones (n, m);
a00=  zeros (n, m);
for i=1:dims
    a0=((source(:,i) * ones(1,m) - ones(n,1) * target(:,i)').^2);
    a000=(source(:,i) * ones(1,m) - ones(n,1) * target(:,i)');
    for j=1:size(Nm,2)
    a00=a00+((source(Nm(:,j),i) * ones(1,m)-a000 - ones(n,1) * target(Nf(:,j),i)').^2);
    end
    M0=M0+a0+size(Nm,2)^2*T*(a00);
    a00=0;
end

if n==m % for non outlier case
    M=round(M0*1e6);
    M=lap(M);
else % for outlier case
    M=round(M0*1e6);
    M=lap_wrapper(M,1e9);
end

end