function [ftc,d,Geod] = GeodesicElasticClosed(p1,p2,stp,figs)

% input p1 and p2 as 2xn or 3xn matrices (planar vs. 3D closed curves)
% to turn off figures set figs=0
% the output is the distance d between p1 and p2 and the geodesic path Geod
% stp refers to the number of shapes displayed along the geodesic

[n,T]=size(p1);


p1 = ReSampleCurve(p1,100);
p2 = ReSampleCurve(p2,100);

q1 = curve_to_q(p1);
q2 = curve_to_q(p2);

tic
[q2n,R] = Find_Rotation_and_Seed_unique(q1,q2,1);
q2n = q2n/sqrt(InnerProd_Q(q2n,q2n));
q2=curve_to_q(p2);
q2 = ProjectC(q2);
p2=q_to_curve(q2);
p2n=R*p2;
toc

d = acos(InnerProd_Q(q1,q2n));
if figs
alpha = geodesic_sphere_Full(q1,q2n,stp);
ftc = Path_Plot(alpha,p2n,10,'b',[73,6],figs);
axis xy;
end

[k] = size(alpha,3);

for j=1:k
    Geod(:,:,j)=q_to_curve(alpha(:,:,j));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
