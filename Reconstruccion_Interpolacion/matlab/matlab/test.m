function ee=test(poly)

[~,n]=size(poly);

ee=zeros(n,2);

for i=1:n
    ee(i,:)=-poly(:,i)+poly(:,mod(i,n)+1);
end
return