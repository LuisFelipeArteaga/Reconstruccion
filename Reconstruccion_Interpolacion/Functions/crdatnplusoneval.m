% % Evaluate cubic Cardinal spline at n+1 values for given four points
% % and tesion. Uniform parameterization is used.

% % INPUT
% % P0,P1,P2,P3:  Given four points in N-dimional space such that
% %               P1 & P2 are endpoints of spline.
% %               P0 & P3 are used to calculate the slope of the endpoints
% %               (i.e slope of P1 & P2).
% %               For example: for 2-dimensional data P0=[x, y],
% %               For example: for 3-dimensional data P0=[x, y, z]
% %               similar analogy for P1, P2, AND P3

% % T: Tension (T=0 for Catmull-Rom type)
% % n: number of intervals (spline is evaluted at n+1 values).

% % OUTPUT
% % MatNbyNPlusOne: Evaluated (interpolated) values of cardinal spline at 
% %                 parameter value u. Values are in N-dimensional space.

function [MatNbyNPlusOne]=crdatnplusoneval(points,T,n)

MatNbyNPlusOne=[];

% u vareis b/w 0 and 1.
% at u=0 cardinal spline reduces to P1.
% at u=1 cardinal spline reduces to P2.

u=0;
MatNbyNPlusOne(:,1)=[evalcrdnd(P0,P1,P2,P3,T,u)]'; % MatNbyNPlusOne(:,1)=length(P0)
du=1/n;
for k=1:n
    u=k*du;
    MatNbyNPlusOne(:,k+1)=[evalcrdnd(P0,P1,P2,P3,T,u)]';
end


function [Pu] =evalcrdnd(P0,P1,P2,P3,T,u)

Pu=[];

s= (1-T)./2;
% MC is cardinal matrix
MC=[-s     2-s   s-2        s;
    2.*s   s-3   3-(2.*s)   -s;
    -s     0     s          0;
    0      1     0          0];

for i=1:length(P0)
    G(:,i)=[P0(i);   P1(i);   P2(i);   P3(i)];
end

U=[u.^3    u.^2    u    1];

for i=1:length(P0)
    Pu(i)=U*MC*G(:,i);
end






