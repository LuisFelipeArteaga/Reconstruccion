function [k,c,angles]=computeCurvature(E)

[n,~]=size(E);

l=zeros(n,1);
angles=zeros(n,1);
k=zeros(n,1);
c=zeros(n,1);
t=zeros(n,2);

for i=1:n
  l(i)=norm(E(i,:));
  t(i,:)=E(i,:)/l(i);
end

for i=1:n
    im=mod(i-2+2*n,n)+1;
    
    cross=abs(det([E(i,:);E(im,:)]));
    dott=dot(E(i,:),E(im,:));
    
    angles(i)=pi-atan2(cross,dott);
    c(i)=2/(l(i)+l(im));
    k(i)=-(angles(i)-pi)*c(i);
end

% for i from 1 to n do
%   im:=(i - 2 + 2*n) mod n +1;
%   ip:=i mod n +1;
% 
%   deltax:=tx[i]-tx[im];
%   deltay:=ty[i]-ty[im];
%   num:=sqrt(deltax^2+deltay^2);
% 
%   deltax:=vx[ip]-vx[im];
%   deltay:=vy[ip]-vy[im];
%   ddenom:=sqrt(deltax^2+deltay^2);
% 
% 
%   k[i]:=2*num/ddenom;
% od:


return