function plotPolygon(poly,showNumbers, cc)


if nargin>2
    col=cc;
else
    col='k';
end


plot([poly(:,1); poly(1,1)],[poly(:,2);poly(1,2)],'o-black',...
    'LineWidth',1,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b',...
    'Color',col,...
    'MarkerSize',4);

axis equal;

if nargin>1 && showNumbers ==1
    [nVertices,~]=size(poly);
   for i=1:nVertices
       text(poly(i,1)+0.01, poly(i,2), num2str(i), 'Color', 'black');
   end
end

return;