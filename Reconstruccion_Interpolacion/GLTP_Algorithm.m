function [Transform, C, T,Y]=GLTP_Algorithm(X, Y, opt)
[N, D]=size(X); [M, D]=size(Y);
% Normalizacion
Yo = Y;
if opt.normalize
    [X,Y,normal] = cpd_normalize(X,Y);
end
[C, W, sigma2, iter, T] =cpd_GRBF(X, Y, opt);

Transform.iter=iter;
Transform.method='cpd';
Transform.Y=T;

Transform.beta=opt.beta;
Transform.W=W;
Transform.Yo=Yo;
Transform.s=1;
Transform.t=zeros(D,1);
if opt.normalize
    Transform.normal=normal;
    Y = T*normal.xscale+repmat(normal.xd,M,1);
    Transform.Y=T*normal.xscale+repmat(normal.xd,M,1);
end

end
%% CPD
function [C, W, sigma2, iter, T] =cpd_GRBF(X, Y, opt)
%Parametros
beta  = opt.beta;
lambda = opt.lambda;
max_it  = opt.max_it;
tol  = opt.tol;

outliers = opt.outliers;

corresp = opt.corresp;
sigma2 = opt.sigma2;
K = opt.K;
alpha = opt.alpha;

[N, D]=size(X); [M, D]=size(Y);
%% Inicicalizacion
% W
T = Y;
W = zeros(M,D);
% Calcular sigma
if ~exist('sigma2','var') || isempty(sigma2) || (sigma2==0)
    sigma2 = (M*trace(X'*X)+N*trace(Y'*Y)-2*sum(X)*sum(Y)')/(M*N*D);
end

%%  Calcular Kernel Gaussiano
G=cpd_G(Y,Y,beta);

%% LLE and Calcular M
d = K;
[~,Lw] = LLE(X',K,d);

MM = (eye(M) - Lw')*(eye(M) -Lw')';
%%  EM Optimizacion
iter  = 0; ntol = tol +10; L = 1;
while (iter<max_it) && (ntol > tol) && (sigma2 > 1e-8)
    
    %% E-step calcular P
    L_old=L;
    [P1,Pt1, PX, L]=cpd_P(X,T, sigma2 ,outliers);
    L=L+lambda/2*trace(W'*G*W);
    ntol=abs((L-L_old)/L);
    disp([' CPD nonrigid  iter= ' num2str(iter), 'sigma '  num2str(sigma2)]);
    %% M-step
    dP = spdiags(P1,0,M,M);
    %     W=(dP*G+lambda*sigma2*eye(M))\(PX-dP*Y);
    %idP=spdiags(1./P1,0,M,M);
    W=(dP*G+sigma2*alpha*eye(M)+sigma2*lambda*MM*G)\(PX-(dP+ sigma2*lambda*MM)*T);
    %Actualizar posiciones
    T=Y+G*W;
    Np=sum(P1);
    %sigma2=abs((sum(sum(X.^2.*repmat(Pt1,1,D)))+sum(sum(T.^2.*repmat(P1,1,D))) -2*trace(PX'*T)) /(Np*D));
    dPt = spdiags(Pt1,0,M,M);
    sigma2 = (trace(X'*dPt*X) -2*trace(T'*PX) - 2*trace(W'*G'*PX) + trace(T'*dP*T) ...
        + 2*trace(W'*G'*dP*T)+trace(W'*G'*dP*G*W))/(Np*D);
    
    pause
    %% Plot iteraciones
    if 1
        figure(1);clf; hold on; axis equal
        scatter(X(:,1),X(:,2),'k')
        scatter(T(:,1),T(:,2),10,'r')
        drawnow
        hold off
    end
    %%
    iter=iter+1;
end

if ~corresp
    C = cpd_corresp(X,T);
end

end
%%
%Normaliza la  conjunto de puntos de referencia and template  a varianza
%media cero y varianza 1
function  [X, Y, normal] =cpd_normalize(x,y)
n=size(x,1);
m=size(y,1);

normal.xd=mean(x);
normal.yd=mean(y);

x=x-repmat(normal.xd,n,1);
y=y-repmat(normal.yd,m,1);

normal.xscale=sqrt(sum(sum(x.^2,2))/n);
normal.yscale=sqrt(sum(sum(y.^2,2))/m);

X=x/normal.xscale;
Y=y/normal.yscale;
end

%%  Gaussian affinity matrix
function G=cpd_G(x,y,beta)
k=-2*beta*beta;
[n, d]=size(x); [m, d]=size(y);

G = pdist2(x,y).^2;
G=exp(G/k);

end
%% Calcular P

function  [P1,Pt1, PX,L] = cpd_P(X,T, sigma2 ,outliers)
[N, D]=size(X); [M, D]=size(T);
ksig = -2*sigma2;
c = (outliers*M*((2*pi*sigma2)^(D/2))) / ((1- outliers)*N);

K = pdist2(X,T).^2;
K = exp(K/ksig);

sK = sum(K,2) +c;
Pt1 = (1- c./sK);
P1 = sum(K./repmat(sK,1,M))' ;
PX = K'*(X./repmat(sK,1,D));
L  = -sum(log(sK));
end

%% Correspondencia entre Y ans X;
function  [C] = cpd_corresp(X,Y)
[N, D]=size(X); [M, D]=size(Y);
if M >= N
    C = zeros(N,1);
    for ii= 1:N
        auxY = Y;
        while true
            [~,posMd] =min(pdist2(X(ii,:),auxY));
            if isempty(find(C == posMd,1))
                C(ii) = posMd;
                break;
            else
                auxY(posMd,:) = NaN;
            end
        end
    end
    
else
    C = zeros(M,1);
    for ii= 1:M
        auxX = X;
        while true
            [~,posMd] =min(pdist2(Y(ii,:),auxX));
            if isempty(find(C == posMd,1))
                C(ii) = posMd;
                break;
            else
                auxX(posMd,:) = NaN;
            end
        end
    end
end

end

