%% Construction of Smooth Branching Surfaces using T-splines

close all; clear; clc

%%  Control contornos iniciales

t = linspace(0,2*pi,100);
c1x = 1.4*cos(t); c1y = 1.1*sin(t);
c2x = .2*cos(t)-.4; c2y = .15*sin(t)+.1;
c3x = .1*cos(t)+.2; c3y = .12*sin(t)-.1;
c4x = .2*cos(t)-.1; c4y = .1*sin(t)+.5;
c5x = .2*cos(t)+.6; c5y = .17*sin(t)-.6;
% Plano 0
Co1 = [c1x;c1y];
% Plano 1
nC = 4; % numero de contornos
Ci1 = [c2x;c2y];mCi1 = mean(Ci1,2);
Ci2 = [c3x;c3y];mCi2 = mean(Ci2,2);
Ci3 = [c4x;c4y];mCi3 = mean(Ci3,2);
Ci4 = [c5x;c5y];mCi4 = mean(Ci4,2);
figure(1);clf; hold on;
plot(Co1(1,:),Co1(2,:),'.-b');
plot(Ci1(1,:),Ci1(2,:),'.-r');
plot(Ci2(1,:),Ci2(2,:),'.-r');
plot(Ci3(1,:),Ci3(2,:),'.-r');
plot(Ci4(1,:),Ci4(2,:),'.-r');

%% Control Cage
%Generacion de Control Network
% Step 0: discretizacion de los contornos iniciales
%% Step 1
% Convex-hull plano 1
piX = [Ci1(1,:) Ci2(1,:) Ci3(1,:) Ci4(1,:) ];
piY = [Ci1(2,:) Ci2(2,:) Ci3(2,:) Ci4(2,:)];
pXY = [piX;piY];
%piX = [Ci1(1,:) Ci2(1,:) Ci3(1,:) ];
%piY = [Ci1(2,:) Ci2(2,:) Ci3(2,:) ];
[k,V] = convhull(piX,piY);
auxK = k(1:end-1);
vCH = [];
pointsC = {};

for ii = 1:nC
    pos = sort(find((auxK/(ii*100)<=1)==1));
    if numel(pos) ~= 0
        pointsC{ii,1} = auxK(pos);
        auxK(pos) = [];
    else
        pointsC{ii,1} = NaN;
    end
    
end
close all;
figure(2);  hold on, axis equal
plot(piX(k),piY(k),'.-r',piX,piY,'+')

%% Step 2
% Voroni triangles identification


%voronoi(piX,piY)
TRI = delaunay(piX,piY);
Tc = zeros(size(TRI,1),nC);
for ii = 1:nC
    if ii == 1
        Tc(:,ii) = sum(TRI>=1 & TRI<= 100,2)==1;
    else
        Tc(:,ii) = sum(TRI>=(ii-1)*100 & TRI<= ii*100,2)==1;
    end
end
posT = find(sum(Tc,2) ==3);

close all;
trimesh(TRI(posT,:),piX,piY); hold on; plot(piX(k),piY(k),'.-r',piX,piY,'.')


%% Step 3
% Modification of contour connectors

% Triangulos vertices
vt = TRI(posT,:);
% vertices triangulos
mA = zeros(numel(vt)/3,3);
for ii = 1:numel(vt)/3
    v1 = [piX(vt(ii,1)),piY(vt(ii,1))];
    v2 = [piX(vt(ii,2)),piY(vt(ii,2))];
    v3 = [piX(vt(ii,3)),piY(vt(ii,3))];
    [aa,ab,ac] = Angulo_T(v1,v2,v3);
    mA(ii,:) =  [aa,ab,ac];
end



%% Reemplazar Convex Hull con los lados internos de Voronoi triangle
close all;
figure(3); hold on
alpha  = deg2rad(120);
[Tri,pVert] = find(mA >= alpha == 1);
nn = [];
vertU = [];
for pp = 1:numel(pVert)
    switch pVert(pp)
        case 1
            ii = Tri(pp); nn = [ii nn];
            v1 = [piX(vt(ii,1)),piY(vt(ii,1))];
            v2 = [piX(vt(ii,2)),piY(vt(ii,2))];
            
            [~,posm1] = min(sqrt(sum(([piX(k)' piY(k)']-v1).^2,2)));
            [~,posm2] = min(sqrt(sum(([piX(k)' piY(k)']-v2).^2,2)));
            if k(posm1)<vt(ii,1)
                newK1 = k(posm1):vt(ii,1);
            else
                newK1 = sort(k(posm1):-1:vt(ii,1));
            end
            if k(posm2)<vt(ii,2)
                newK2 = k(posm2):vt(ii,2);
            else
                newK2 = sort(k(posm2):-1:vt(ii,2));
            end
            vertU = [[piX(vt(ii,3)),piY(vt(ii,3))]; vertU];
            k = [k(1:posm1-1); newK1'; vt(ii,3); newK2'; k(posm1+1:end) ];
        case 2
            ii = Tri(pp); nn = [ii nn];
            v2 = [piX(vt(ii,2)),piY(vt(ii,2))];
            v3 = [piX(vt(ii,3)),piY(vt(ii,3))];
            
            [~,posm1] = min(sqrt(sum(([piX(k)' piY(k)']-v2).^2,2)));
            [~,posm2] = min(sqrt(sum(([piX(k)' piY(k)']-v3).^2,2)));
            if k(posm1)<vt(ii,2)
                newK1 = k(posm1):vt(ii,2);
            else
                newK1 = sort(k(posm1):-1:vt(ii,2));
            end
            if k(posm2)<vt(ii,3)
                newK2 = k(posm2):vt(ii,3);
            else
                newK2 = sort(k(posm2):-1:vt(ii,3));
            end
            vertU = [[piX(vt(ii,1)),piY(vt(ii,1))];vertU];
            k = [k(1:posm1-1); newK1'; vt(ii,1); newK2'; k(posm1+1:end) ];
        case 3
            ii = Tri(pp);  nn = [ii nn];
            v1 = [piX(vt(ii,1)),piY(vt(ii,1))];
            v3 = [piX(vt(ii,3)),piY(vt(ii,3))];
            
            [~,posm1] = min(sqrt(sum(([piX(k)' piY(k)']-v1).^2,2)));
            [~,posm2] = min(sqrt(sum(([piX(k)' piY(k)']-v3).^2,2)));
            if k(posm1)<vt(ii,1)
                newK1 = k(posm1):vt(ii,1);
            else
                newK1 = sort(k(posm1):-1:vt(ii,1));
            end
            if k(posm2)<vt(ii,3)
                newK2 = k(posm2):vt(ii,3);
            else
                newK2 = sort(k(posm2):-1:vt(ii,3));
            end
            vertU = [[piX(vt(ii,2)),piY(vt(ii,2))];vertU];
            k = [k(1:posm1-1); newK1'; vt(ii,2); newK2'; k(posm1+1:end) ];
    end
end
vt(nn,:) = [];
% trimesh(TRI(posT(1),:),piX,piY);
% plot(piX(k),piY(k),'.-r',piX,piY,'.')
% plot(mCi1(1),mCi1(2),'+');text(mCi1(1),mCi1(2),'C1')
% plot(mCi2(1),mCi2(2),'+');text(mCi2(1),mCi2(2),'C2')
% plot(mCi3(1),mCi3(2),'+');text(mCi3(1),mCi3(2),'C3')
% plot(mCi4(1),mCi4(2),'+');text(mCi4(1),mCi4(2),'C4')


%% Step 4
% connecting Voronoi vertices to contours
% Baricentro

for qq =1: numel(vt)/3
    v1 = [piX(vt(qq,1)),piY(vt(qq,1))];
    v2 = [piX(vt(qq,2)),piY(vt(qq,2))];
    v3 = [piX(vt(qq,3)),piY(vt(qq,3))];
    bcT =[(v1(1)+v2(1)+v3(1))/3 (v1(2)+v2(2)+v3(2))/3];
    
    
    plot(bcT(1),bcT(2),'*')
    plot([bcT(1) v1(1)],[bcT(2) v1(2)],'.-b')
    plot([bcT(1) v2(1)],[bcT(2) v2(2)],'.-b')
    plot([bcT(1) v3(1)],[bcT(2) v3(2)],'.-b')
end


triC=ceil(vt/100);
vCHi = [1; find(diff(ceil(k/100))); numel(k)];
ki = (numel(k)+1) - find(diff(ceil(flipud(k)/100)));
vCH = unique([ki;vCHi]);
kCH = k(vCH(1:end-1));

vertC = ceil(kCH/100);
%vertC = [vertC(1:4);5; 5; vertC(5:end) ];
isCH = find(ismember(vertC,triC) == 0);
isCH = [isCH(1)-1 ; isCH; isCH(end)+1 ];
D = diff([0,diff(isCH')==1,0]);
Di = find((D>0)); Df = find((D<0));
vertL = [];
%% struct linea entre contornos
lineCont = struct('p1',[nan nan],'pmean',[nan nan],'p2',[nan nan],'cont', [nan nan]);
%%
posMM= [];
for qq = 1:numel(D)/2-1
    auxD = kCH(isCH(Di(qq):Df(qq)));
    newV = [mean(piX(auxD)),mean(piY(auxD))];
    conV = unique(vertC(isCH(Di(qq):Df(qq))));
    
    lineCont(qq).pmean = newV;
    for pp = 1:numel(conV)
        auxC = [piX(((conV(pp)-1)*100)+1:conV(pp)*100);piY(((conV(pp)-1)*100)+1:conV(pp)*100)];
        [~,posM]= min(sqrt(sum((auxC-newV').^2)));
        posM = ((conV(pp)-1)*100)+1+posM;
        posMM= [posMM; posM];
        plot([ newV(1) piX(posM)],[newV(2) piY(posM)],'*-b')
        vertL = [vertL;[piX(posM) piY(posM)];newV]
    end
    lineCont(qq).p1 = [piX(posMM(1)) piY(posMM(1))];
    lineCont(qq).p2 = [piX(posMM(2)) piY(posMM(2))];
    lineCont(qq).cont =  ceil(posMM/100)';
    plot(newV(1),newV(2),'*')
end


% Separra vertice comun
vertices = [piX(kCH);piY(kCH)];
a = zeros(2,size(vertU,2));
radV = 0.04;
aCell = cell(2,size(vertU,2));
aux = 1;
for qq = 1:size(vertU,2)
    posV = find( sum([piX(kCH);piY(kCH)] == vertU(qq,:)')==2);
    pp = vertC(posV);
    plot(vertices(1,posV),vertices(2,posV),'*m')
    auxC = [piX(((pp-1)*100)+1:pp*100);piY(((pp-1)*100)+1:pp*100)];
    distV = sqrt(sum((auxC-vertices(:,posV)).^2));
    a = find((distV<radV)) + ((pp-1)*100)+1;
    aCell{qq,1} = a;
    aCell{qq,2} = posV;
end

nk =[];

auxvCH =[0 ;vCH([aCell{:,2}]);numel(k)+1];
for qq = 1:size(aCell,1)+1
    if qq==size(aCell,1)+1
        auxK = [k(auxvCH(qq)+1:auxvCH(qq+1)-1)];
        nk = [nk;auxK];
    else
        auxK = [k(auxvCH(qq)+1:auxvCH(qq+1)-1);aCell{qq,1}'];
        nk = [nk;auxK];
    end
end

vCHi = [1; find(diff(ceil(nk/100))); numel(nk)];
ki = (numel(nk)+1) - find(diff(ceil(flipud(nk)/100)));
vCH = unique([ki;vCHi]);
kCH = nk(vCH(1:end-1));

plot(piX(kCH),piY(kCH),'*')
trimesh(TRI(posT(1),:),piX,piY);
plot(piX(nk),piY(nk),'.-r',piX,piY,'.')
plot(mCi1(1),mCi1(2),'+');text(mCi1(1),mCi1(2),'C1')
plot(mCi2(1),mCi2(2),'+');text(mCi2(1),mCi2(2),'C2')
plot(mCi3(1),mCi3(2),'+');text(mCi3(1),mCi3(2),'C3')
plot(mCi4(1),mCi4(2),'+');text(mCi4(1),mCi4(2),'C4')

%%  Step 5- surrounding curve(SC)m,
close all

poligLine = [kCH,[kCH(2:end);kCH(1)]];

LineI = pXY(:,poligLine(:,1));
LineF = pXY(:,poligLine(:,2));

figure(2);clf
hold on , axis equal
for pp = 1:numel(poligLine)/2
    %     plot([LineI(1,pp) LineF(1,pp) ] ,[LineI(2,pp) LineF(2,pp) ],'o-k')
end
plot(piX,piY,'.')
plot(mCi1(1),mCi1(2),'+');text(mCi1(1),mCi1(2),'C1')
plot(mCi2(1),mCi2(2),'+');text(mCi2(1),mCi2(2),'C2')
plot(mCi3(1),mCi3(2),'+');text(mCi3(1),mCi3(2),'C3')
plot(mCi4(1),mCi4(2),'+');text(mCi4(1),mCi4(2),'C4')
%Plano 0
plot(c1x,c1y,'-b','LineWidth',2)


%
vertSC = pXY(:,kCH)';
Co1 = Co1';

%[xy,distance,t_a]= distance2curve(vertSC,Co1,'linear');
hold on
% plot(Co1(:,1),Co1(:,2),'k.')
%% median vertices Step 5

%%%%%%%%%%%%%%%% Struct Vertices
pointsV = struct('xy',[nan nan],'type','type','cont',nan,'ind',nan,'C',nan,'Lpm',[nan nan]);
kCHn = ceil(kCH/100);
for ii =1:numel(kCH)
    pointsV(ii).xy = vertSC(ii,:);
    pointsV(ii).type = 'VerV';
    pointsV(ii).cont = kCHn(ii);
    pointsV(ii).ind =ii;
end
%%%%%%%%%%%%%%%%
pmC = find(diff (ceil(kCH/100)));
pm = zeros(numel(pmC),2);
for qq = 1:numel(pmC)
    p1 = vertSC(pmC(qq),:);
    p2 = vertSC(pmC(qq)+1,:);
    pm(qq,:) = p1*.5+(1-0.5)*p2;
    %      plot(pm(qq,1),pm(qq,2),'.r')
end
qq  = qq +1;
p1 = vertSC(1,:);
p2 = vertSC(end,:);
pm(qq,:) = p1*.5 + .5*p2;
% plot(pm(qq,1),pm(qq,2),'.r')

vertSCN = [];
npmC = [0; pmC ; numel(kCH)];
for ii = 1:numel(npmC)-1
    auxV = vertSC(npmC(ii)+1:npmC(ii+1),:);
    vertSCN = [vertSCN;auxV;pm(ii,:)];
end
pmC = npmC(2:end)
%%%%%%%%%%%%%%%%%%% Struct PM
pointsPM = struct('xy',[nan nan],'type','type','betw',[nan nan],'cont',[nan nan],'ind',nan,'C',nan);
kCHn = ceil(kCH/100);
for ii =1:numel(pm)/2
    pointsPM(ii).xy = pm(ii,:);
    pointsPM(ii).type = 'pm';
    if pmC(ii) + 1 > pmC(end)
        pointsPM(ii).betw= [pmC(ii) 1];
        pointsPM(ii).cont = [kCHn(pmC(ii)) kCHn(1)];
    else
        pointsPM(ii).betw= [pmC(ii) pmC(ii)+1];
        pointsPM(ii).cont = [kCHn(pmC(ii)) kCHn(pmC(ii)+1)];
    end
    pointsPM(ii).ind =ii;
end
%% CPD - Linking   SC   to   contour   on plane-0
opt.method = 'nonrigid';opt.normalize = 0;
opt.max_it = 100; opt.tol = 1e-5;
opt.viz = 1; opt.corresp = 1;
opt.outliers = .01; opt.fgt = 0;
opt.sigma2 = 1;
% reflejan la cantidad de regularización de la suavidad.
opt.beta = .1; opt.lambda = 20;

vertSCN = [vertSCN; vertSCN(1,:)];
smoothV = interparc(linspace(0,1,100),vertSCN(:,1),vertSCN(:,2));
smoothV=smoothV(1:end-1,:);


[~,bm]=min(pdist2(vertSCN(1:end-1,:),smoothV),[],2);

addpath('mex/')
[Transform, C]=cpd_register(Co1, smoothV, opt);
smoothV = smoothV(bm,:);
C = C(bm,:);

Cpm =  C(3:3:end)
for ii = 1:numel(Cpm)
    pointsPM(ii).C= Cpm(ii);
end
C1 = C;
C1(3:3:end)= [];
for ii = 1: numel(C1)
    pointsV(ii).C = C1(ii);
end
%%
clc
figure(2); hold on

% for pp = 1:5
%     plot([pointsPM(pp).xy(1) Co1(pointsPM(pp).C,1) ] ,...
%         [pointsPM(pp).xy(2) Co1(pointsPM(pp).C,2)],'b')
% end
for pp = 1:10
    plot([pointsV(pp).xy(1) Co1(pointsV(pp).C,1) ] ,...
        [pointsV(pp).xy(2) Co1(pointsV(pp).C,2)],'k','LineWidth',2)
end
% plot vertices triangulo
%%
trianV = struct('v1' ,[ nan nan],'v2' , [nan nan],'v3', [nan nan],'cent' ,[nan nan],'cont', [1 2 3]);

%%

for qq =1: numel(vt)/3
    v1 = [piX(vt(qq,1)),piY(vt(qq,1))];
    v2 = [piX(vt(qq,2)),piY(vt(qq,2))];
    v3 = [piX(vt(qq,3)),piY(vt(qq,3))];
    trianV(qq).v1  = v1; trianV(qq).v2  = v2; trianV(qq).v3  = v3;
    trianV(qq).cont = ceil(vt/100);
    bcT =[(v1(1)+v2(1)+v3(1))/3 (v1(2)+v2(2)+v3(2))/3];
    trianV(qq).cent  = bcT;
    
    %plot(bcT(1),bcT(2),'*')
    plot([bcT(1) v1(1)],[bcT(2) v1(2)],'k','LineWidth',2)
    plot([bcT(1) v2(1)],[bcT(2) v2(2)],'k','LineWidth',2)
    plot([bcT(1) v3(1)],[bcT(2) v3(2)],'k','LineWidth',2)
end
%
for ii = 1:1
    plot([lineCont(ii).p1(1) lineCont(ii).pmean(1) lineCont(ii).p2(1)],...
        [lineCont(ii).p1(2) lineCont(ii).pmean(2) lineCont(ii).p2(2)],'-k','LineWidth',2 )
end

%% Linking Voronoi vertices

clearvars -except trianV pointsV pointsPM lineCont Co1 Ci1 Ci2 Ci3 Ci4
% ------------ Triangulo----------------
for ii = 1:numel(trianV)
    pPM = reshape([pointsPM(:).cont],[2,numel(pointsPM)])';
    vV = nchoosek(trianV(ii).cont,numel(trianV(ii).cont)-1);
    for pp = 1:3
        pnV = find(sum(ismember(pPM,vV(pp,:)),2) ==2);
        pointsPM(pnV).xy = trianV.cent;
    end
end

% -------------Line-----------------

for ii = 1:numel(trianV)
    pPM = reshape([pointsPM(:).cont],[2,numel(pointsPM)])';
    vV =lineCont(ii).cont;
    
    pnV = find(sum(ismember(pPM,vV),2) ==2);
    for pp = 1:numel(pnV)
        pointsPM(pnV(pp)).xy = lineCont.pmean;
    end
end


for pp = 1:5
    plot([pointsPM(pp).xy(1) Co1(pointsPM(pp).C,1) ] ,...
        [pointsPM(pp).xy(2) Co1(pointsPM(pp).C,2)],'r','LineWidth',2)
end
%% Step 7 bridge-edges

%----------------------- 1-1

% ----------------------- 0-1a
pointsV(:).cont;

for ii = 1:2:numel(pointsV)-1
    Lpm1 = .5*pointsV(ii).xy + .5*[Co1(pointsV(ii).C,1) Co1(pointsV(ii).C,2)];
    Lpm2 = .5*pointsV(ii+1).xy + .5*[Co1(pointsV(ii+1).C,1) Co1(pointsV(ii+1).C,2)];
    pointsV(ii).Lpm = Lpm1;
    pointsV(ii+1).Lpm = Lpm2;
    plot([Lpm1(1) Lpm2(1)],[Lpm1(2) Lpm2(2)],'k','LineWidth',2)
end
% ----------------------- 0-1b
% Linking triangulo
for ii =1:numel(trianV)
    pPM = reshape([pointsPM(:).cont],[2,numel(pointsPM)])';
    pnV = find(sum(ismember(pPM,trianV(ii).cont),2) ==2);
    LinkC1 = [pointsPM(pnV).betw];
    for pp = 1:numel(LinkC1)
        plot([pointsV(LinkC1(pp)).Lpm(1)  trianV(ii).cent(1)],...
            [pointsV(LinkC1(pp)).Lpm(2)  trianV(ii).cent(2)],'k','LineWidth',2)
    end
end
%linking Triangulo
for ii =1:numel(lineCont)
    pPM = reshape([pointsPM(:).cont],[2,numel(pointsPM)])';
    pnV = find(sum(ismember(pPM,lineCont(ii).cont),2) ==2);
    LinkC2 = [pointsPM(pnV).betw];
    for pp = 1:numel(LinkC2)
        plot([pointsV(LinkC2(pp)).Lpm(1)  lineCont(ii).pmean(1)],...
            [pointsV(LinkC2(pp)).Lpm(2)  lineCont(ii).pmean(2)],'k','LineWidth',2)
    end
end

%% plot Contours
%Plano 0
Ct = sort([pointsV(:).C pointsPM(:).C]); Ct = [Ct Ct(1)];
plot(Co1(Ct,1),Co1(Ct,2),'--b','LineWidth',1.5)

% Contornos plano 1
sPlano1= struct('cont',nan,'xy' ,[nan nan]);
cV = [pointsV.cont];
cPM = [trianV.cont];
cLC = [lineCont.cont];
for ii = 1:4
    vcPM = [];  vcLC = [];
    fcV = find(cV==ii);
    vcV =cell2mat({pointsV(fcV).xy}');
    fcPM = find(cPM==ii);
    if fcPM == 1;vcPM = (trianV.v1);elseif  fcPM == 2;vcPM = (trianV.v2);elseif  fcPM == 3;vcPM = (trianV.v3);end
    fcLC = find(cLC==ii);
    if fcLC== 1;vcLC = (lineCont.p1);elseif fcLC== 2;vcLC = (lineCont.p2);end
    cVT = [vcV;vcPM;vcLC];
    K = convhull(cVT(:,1),cVT(:,2));
    plot(cVT(K,1),cVT(K,2),'-.b','LineWidth',1.5)
    
    sPlano1(ii).cont = ii;
    sPlano1(ii).xy = cVT(K,:);
end

%% ----------------------------------Cortornos 3D 
clearvars -except trianV pointsV pointsPM lineCont Co1 Ci1 Ci2 Ci3 Ci4  LinkC2  LinkC1 sPlano1
figure(3); axis equal; clf; hold on;grid on
%axis([-2 2 -2 2 0 1])
clc
z0 = 0;
%---------------------------------Plano 0
Ct = sort([pointsV(:).C pointsPM(:).C]); Ct = [Ct Ct(1)]; 
Co1 = [Co1,z0*ones(size(Co1,1),1)];
plot3(Co1(Ct,1),Co1(Ct,2),Co1(Ct,3),'-.b','LineWidth',1.5)

%---------------------------------Plano 1
z1 =.4;
for  ii= 1:numel(sPlano1)
    sPlano1(ii).xy = [sPlano1(ii).xy, z1*ones(size(sPlano1(ii).xy,1),1)];
    plot3(sPlano1(ii).xy(:,1),sPlano1(ii).xy(:,2),sPlano1(ii).xy(:,3),'-.r','LineWidth',1.5)
end

%--------------------------------Plano medio
zm =.5*(z0+z1);
% Linking vertices-Plano 0
for ii = 1:numel(pointsV)
    pointsV(ii).xy = [pointsV(ii).xy z1];
    plot3([pointsV(ii).xy(1) Co1(pointsV(ii).C,1)],[pointsV(ii).xy(2) Co1(pointsV(ii).C,2)],...
        [pointsV(ii).xy(3) Co1(pointsV(ii).C,3)],'k','LineWidth',2)
end
% Linking puntos medios
for ii = 1:2:numel(pointsV)-1
    pointsV(ii).Lpm = [pointsV(ii).Lpm zm];
    pointsV(ii+1).Lpm = [pointsV(ii+1).Lpm zm];
    plot3([pointsV(ii).Lpm(1) pointsV(ii+1).Lpm(1)],[pointsV(ii).Lpm(2) pointsV(ii+1).Lpm(2)],...
        [pointsV(ii).Lpm(3) pointsV(ii+1).Lpm(3)],'k','LineWidth',2)
end

% Linking centroide-plano 0 
for pp = 1:5
   pointsPM(pp).xy = [pointsPM(pp).xy zm] ;
    plot3([pointsPM(pp).xy(1) Co1(pointsPM(pp).C,1) ] ,[pointsPM(pp).xy(2) Co1(pointsPM(pp).C,2)],...
        [pointsPM(pp).xy(3) Co1(pointsPM(pp).C,3) ],'k','LineWidth',2)
end

% Linking centroide- plano1
for ii = 1:numel(trianV)
    trianV(ii).v1 = [trianV.v1 z1];
    trianV(ii).v2 = [trianV.v2 z1];
    trianV(ii).v3 = [trianV.v3 z1];
    trianV(ii).cent = [trianV(ii).cent zm];
    
    
    plot3([trianV(ii).v1(1) trianV(ii).cent(1)],[trianV(ii).v1(2) trianV(ii).cent(2)],[trianV(ii).v1(3) trianV(ii).cent(3)],'k','LineWidth',2)
    plot3([trianV(ii).v2(1) trianV(ii).cent(1)],[trianV(ii).v2(2) trianV(ii).cent(2)],[trianV(ii).v2(3) trianV(ii).cent(3)],'k','LineWidth',2)
    plot3([trianV(ii).v3(1) trianV(ii).cent(1)],[trianV(ii).v3(2) trianV(ii).cent(2)],[trianV(ii).v3(3) trianV(ii).cent(3)],'k','LineWidth',2)
end

% Linking line- plano1
for ii = 1:numel(lineCont)
    lineCont(ii).p1 = [lineCont(ii).p1 z1]; 
    lineCont(ii).p2 = [lineCont(ii).p2 z1] ;
    lineCont(ii).pmean = [lineCont(ii).pmean zm];
    plot3([lineCont(ii).p1(1) lineCont(ii).pmean(1)],[lineCont(ii).p1(2) lineCont(ii).pmean(2)],[lineCont(ii).p1(3) lineCont(ii).pmean(3)],'k','LineWidth',2)
    plot3([lineCont(ii).p2(1) lineCont(ii).pmean(1)],[lineCont(ii).p2(2) lineCont(ii).pmean(2)],[lineCont(ii).p2(3) lineCont(ii).pmean(3)],'k','LineWidth',2)
end

%linking puntos medios-centroides
for pp=1:numel(LinkC1)
    plot3([pointsV(LinkC1(pp)).Lpm(1)  trianV.cent(1)],[pointsV(LinkC1(pp)).Lpm(2)  trianV.cent(2)],...
        [pointsV(LinkC1(pp)).Lpm(3)  trianV.cent(3)],'m','LineWidth',2)
end

%linking puntos medios-(media-linea)
 for pp = 1:numel(LinkC2)
        plot3([pointsV(LinkC2(pp)).Lpm(1)  lineCont.pmean(1)],[pointsV(LinkC2(pp)).Lpm(2)  lineCont.pmean(2)],...
            [pointsV(LinkC2(pp)).Lpm(3)  lineCont.pmean(3)],'r','LineWidth',2)
 end
    
 %% 
 clear
 close
 clc
 load c3.mat
 addpath('Concave Hull')
a2 = [.7785 0.2087;0.3859 -0.3293;0.08178 -.6593;-.08063 -.6283;-.02741 0.1793;.6901 .3271;.7785 0.2087];
a1 = [0.02741 0.1796;
    -0.9854 0.106;
    -0.8467 0.2126;
    0.02741 0.1796];
a4= [0.3859 -0.3293;
           .4511 -0.8272;
           1.008 -0.4352;
           0.3859 -0.3293];
figure(1)
axis equal
 hold on;scatter(Ci3(1,:),Ci3(2,:)) 
 scatter(a2(:,1),a2(:,2),'+')
 scatter(a3(:,1),a3(:,2),'*')
 ain = [a3(end,:);a3];
 pt1 = interparc(100,ain(:,1),ain(:,2),'csape');
 pt2 = interparc(100,a2(:,1),a2(:,2),'csape');
 pt3 = interparc(100,a1(:,1),a1(:,2),'csape');
pt4 = interparc(100,a4(:,1),a4(:,2),'csape'); 
 
 
 scatter(Co1(:,1),Co1(:,2))
 scatter(pt1(:,1),pt1(:,2))
 scatter(pt2(:,1),pt2(:,2))
 scatter(pt3(:,1),pt3(:,2))
 scatter(pt4(:,1),pt4(:,2))

dataset = [pt1(1:end-1,:);pt2(1:end-1,:);pt3(1:end-1,:);pt4(1:end-1,:)];
k = boundary(dataset(:,1),dataset(:,2));
hh = concaveHull(dataset,4);
% % DT = delaunayTriangulation(dataset(:,1),dataset(:,2))
% % IC = incenter(DT);
% % triplot(DT)
% % hold on
% % plot(IC(:,1),IC(:,2),'*r')
% hold on;
plot(hh(:,1),hh(:,2));

bb = interparc(100,hh(:,1),hh(:,2),'linear');
c1 = {Ci3',pt1};
c2 = {Ci2',pt2};
c3 = {Ci1',pt3};
c4 = {Ci4',pt4};

c0 = {bb,Co1(:,[1 2])};


 