function poly=buildPolygonFromIntrinsicDefinition(edges, aa, p0)
close all;
angles=cumsum(pi-aa);

assert(length(edges)==length(angles));
nVerices=length(edges);


l=zeros(nVerices,1);

for i=1:nVerices
    l(i)=norm(edges(i,:));
end

poly=zeros(2,nVerices);

poly(:,1)=p0;
poly(:,2)=p0;
poly(1,2)=poly(1,2)+l(1);

for i=2:nVerices-1
    poly(1,i+1)=poly(1,i)+l(i-1)*cos(angles(i-1));
    poly(2,i+1)=poly(2,i)+l(i-1)*sin(angles(i-1));
end

% plotPolygon(poly',1);


return