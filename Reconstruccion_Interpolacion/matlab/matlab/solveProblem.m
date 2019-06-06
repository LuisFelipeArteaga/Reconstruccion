function [angles,poly]=solveProblem(E,K,pp)
if nargin > 2
    p=pp;
else
    p=1;
end

global lengths;
global edges;
edges=E;
nVerices=length(E);
lengths=zeros(nVerices,1);
for i=1:nVerices
    lengths(i)=norm(E(i,:));
end

close all;

power=pp;

if power==Inf
    power=1;
end

[oldK,C,oldAngles]=computeCurvature(E);

[n,~]=size(K);
angles0=rand(nVerices,1)*2;%zeros(n,1);

poly0=buildPolygonFromIntrinsicDefinition(E,angles0,[0,0]);
poly0=realignPoly(poly0);


problem={};
problem.objective=@(x)(1/2*norm(-(x-pi).*C-K,p)^power);
problem.x0=angles0;
problem.Aineq=[];
problem.bineq=[];
problem.Aeq=ones(1,n);
problem.beq=(n-2)*pi;

problem.solver='fmincon';
problem.options=optimoptions('fmincon');
problem.options.MaxFunEvals = 100000;
% problem.options.Algorithm='sqp';

problem.nonlcon=@nonlinear;



angles = fmincon(problem);

disp(['Error ' num2str(norm(-(angles-pi).*C-K,p))]);
disp(['Old Error ' num2str(norm(oldK-K,p))]);
disp(['Angle avg diff ' num2str(mean(abs(oldAngles-angles)))]);


poly=buildPolygonFromIntrinsicDefinition(E,angles,[0,0]);
poly=realignPoly(poly);

oldPoly=buildPolygonFromIntrinsicDefinition(E,oldAngles,[0,0]);
oldPoly=realignPoly(oldPoly);

plotPolygon(poly',0,'r');
hold on;
plotPolygon(oldPoly',0,'g');
plotPolygon(poly0',0,'b');
legend('new poly','old poly','starting poly');

return


function [c,ceq]=nonlinear(x)
global lengths;
global edges;

c=0;
poly=buildPolygonFromIntrinsicDefinition(edges, x, [0,0]);
ceq=norm(poly(:,end))-lengths(end-1);
return


function res=realignPoly(poly)
centre=mean(poly,2);
res=bsxfun(@minus,poly,centre);

cc=res*res';
[r,~]=eig(cc);
res=r*res;


return