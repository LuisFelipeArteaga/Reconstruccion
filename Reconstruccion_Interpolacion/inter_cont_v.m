%%  interpolation varios
clear ; close all; clc
addpath('Functions','Files','ClosedCurves2D3D','Mesh_voxelisation');
%%
t = linspace(0,2*pi,100);
x = cos(t);
y = sin(t);

c1  = [x;y;0*ones(size(t))] +[2;2;0];
c2 = [.35*x;.25*y;1*ones(size(t))]+[.4 ;.4 ;0] +[2;2;0];
c3 = [.25*x;.35*y;1*ones(size(t))] - [.4 ;.4 ;0] +[2;2;0];
c4 = [.25*x;.35*y;1*ones(size(t))] + [.4 ;-.4 ;0] +[2;2;0];

figure(1); hold on; axis equal
plot3(c1(1,:),c1(2,:),c1(3,:),'.')
plot3(c2(1,:),c2(2,:),c2(3,:),'.')
plot3(c3(1,:),c3(2,:),c3(3,:),'.')
plot3(c4(1,:),c4(2,:),c4(3,:),'.')

%%
std = 10;
[ftc1,~,~,Ro1] = GeodesicElasticClosed2(c1(1:2,:),c2(1:2,:),std-1);
[ftc2,~,~,Ro2] = GeodesicElasticClosed2(c1(1:2,:),c3(1:2,:),std-1);
[ftc3,~,~,Ro3] = GeodesicElasticClosed2(c1(1:2,:),c4(1:2,:),std-1);


%%
t2 = linspace(0,1,std);
mc1 = mean(c1,2);
mc2 = mean(c2,2);
mc3 = mean(c3,2);
mc4 = mean(c4,2);

s1 = (1-repmat(t2,3,1)).*mc1 +  (repmat(t2,3,1)).*mc2;
s2 = (1-repmat(t2,3,1)).*mc1 +  (repmat(t2,3,1)).*mc3;
s3 = (1-repmat(t2,3,1)).*mc1 +  (repmat(t2,3,1)).*mc4;

plot3(s1(1,:),s1(2,:),s1(3,:),'+')
plot3(s2(1,:),s2(2,:),s2(3,:),'+')
plot3(s3(1,:),s3(2,:),s3(3,:),'+')

figure(1); hold on, axis equal

di1 = deltaS(c1(1:2,:)',(ftc1{1}*Ro1));
df1 = deltaS(c2(1:2,:)',(ftc1{end}*Ro1));
di2 = deltaS(c1(1:2,:)',(ftc2{1}*Ro2));
df2 = deltaS(c3(1:2,:)',(ftc2{end}*Ro2));
di3 = deltaS(c1(1:2,:)',(ftc3{1}*Ro3));
df3 = deltaS(c4(1:2,:)',(ftc3{end}*Ro3));

ds1 = (1-repmat(t2,2,1)).*di1' +  (repmat(t2,2,1)).*df1';
ds2 = (1-repmat(t2,2,1)).*di2' +  (repmat(t2,2,1)).*df2';
ds3 = (1-repmat(t2,2,1)).*di3' +  (repmat(t2,2,1)).*df3';


%%
pC = [];
for  ii =1:std
    %figure(1); hold on
    aa1 = (ftc1{ii}*Ro1).*ds1(:,ii)';
    aa2 = (ftc2{ii}*Ro2).*ds2(:,ii)';
    aa3 = (ftc3{ii}*Ro3).*ds3(:,ii)';
    
    Cf1 = horzcat(aa1 , s1(3,ii)*ones(numel(aa1)/2,1));
    Cf1 =Translation_matrix(Cf1',s1(:,ii)');
    %plot3(Cf1(1,:),Cf1(2,:),Cf1(3,:),'.-')
    
    Cf2 = horzcat(aa2 , s1(3,ii)*ones(numel(aa2)/2,1));
    Cf2 =Translation_matrix(Cf2',s2(:,ii)');
   % plot3(Cf2(1,:),Cf2(2,:),Cf2(3,:),'.-')
    
    Cf3 = horzcat(aa3 , s1(3,ii)*ones(numel(aa3)/2,1));
    Cf3 =Translation_matrix(Cf3',s3(:,ii)');
    
    mn = [max(Cf1,[],2), max(Cf2,[],2)]; mn = round(mn(1:2,:)+1);
    m = max(mn(1,:)); n = max(mn(2,:));
    
    b1 = poly2mask(Cf1(1,:)*100,Cf1(2,:)*100,100*m,100*n);
    b2 = poly2mask(Cf2(1,:)*100,Cf2(2,:)*100,100*m,100*n);
    b3 = poly2mask(Cf3(1,:)*100,Cf3(2,:)*100,100*m,100*n);
    overL = sum(sum((b1 & b2) & b3));
    
    figure(2); hold on; axis equal
    if overL > 0
%         aaa =[Cf1(1,:),Cf2(1,:);Cf1(2,:),Cf2(2,:)];
%         aaaz = [Cf1(3,:),Cf2(3,:)];
%         %k= convhull(aaa(1,:),aaa(2,:));
%         k= boundary(aaa(1,:)',aaa(2,:)');
%         plot3(aaa(1,k),aaa(2,k),aaaz(k),'.-')

        bw = (b1 + b2) +b3;
        bw2 = bwperim(bw);
        [row,col] = find(bw2 == 1);
        row = row./100;
        col = col./100;
        zz = Cf1(3,1)*ones(size(row));
        plot3(row,col,zz,'.')
        pC  = [[row col zz];pC];
    else
        bw1 = bwperim(b1);
        [row1,col1] = find(bw1 == 1);
        row1 = row1./100;
        col1 = col1./100;
        zz1= Cf1(3,1)*ones(size(row1));
        plot3(row1,col1,zz1,'.')
         pC  = [[row1 col1 zz1];pC];
        
        
        bw2 = bwperim(b2);
        [row2,col2] = find(bw2 == 1);
        row2 = row2./100;
        col2 = col2./100;
        zz2 = Cf2(3,1)*ones(size(row2));
        plot3(row2,col2,zz2,'.')
         pC  = [[row2 col2 zz2];pC];
        
        bw3 = bwperim(b3);
        [row3,col3] = find(bw3 == 1);
        row3 = row3./100;
        col3 = col3./100;
        zz3= Cf3(3,1)*ones(size(row3));
        plot3(row3,col3,zz3,'.')
         pC  = [[row3 col3 zz3];pC];
        %plot3(Cf1(1,:),Cf1(2,:),Cf1(3,:),'.-')
        %plot3(Cf2(1,:),Cf2(2,:),Cf2(3,:),'.-')
    end    
end




















