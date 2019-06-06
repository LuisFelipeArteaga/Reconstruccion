function [param, model] = cpd(config, model, scene,ctrl_pts)
tic
%%
%%
% figure; hold on
% mm = mean(model);
% ms = mean(scene);
% scatter(model(:,1)-mm(1),model(:,2)-mm(2))
% scatter(scene(:,1)-ms(1),scene(:,2)-ms(2))
%%

sigma = config.init_sigma;
anneal_rate = config.anneal_rate;
outliers =0;
lambda = config.lambda;
max_iter = config.max_iter;
max_em_iter = config.max_em_iter;

tol = config.tol;
EMtol = config.emtol;
[n,d] = size(ctrl_pts);
[m,d] = size(model);

[model, ~,~] = cpd_normalize(model);
[ctrl_pts, ~, ~] = cpd_normalize(ctrl_pts);
[scene, centroid, scale] = cpd_normalize(scene);
model0 = model;

beta = config.beta;
basis = cpd_G(model,ctrl_pts,beta);
kernel = cpd_G(ctrl_pts,ctrl_pts,beta);
%A = inv(basis'*basis+lambda*sigma*sigma*kernel)*basis';
A = basis'*basis+lambda*sigma*sigma*kernel;


param = config.init_param;  % it should always be of size n*d
model = model + basis*param;


%it_total = 1;
%flag_stop = 0;

iter=0; E=1; ntol=tol+10;

%%
figure; hold on; axis equal
scatter(model(:,1),model(:,2))
scatter(scene(:,1),scene(:,2))
%%
while (iter < max_iter) && (ntol > tol)
    EMiter=0; EMtol=tol+10;  % repeat at each termperature.
    model_old = model;
    while (EMiter < max_em_iter) && (EMtol > tol)
        %disp(sprintf('EMiter=%d',EMiter));
        %disp(sprintf('E=%f',E));
        E_old = E;
        % E-step: Given model&scene, update posterior probability matrix P.
        [P,Eu] = cpd_P(model, scene, sigma, outliers);
        % M-step: Given correspondence, solve warp parameter.
        %
        E = Eu + lambda/2*trace(param'*kernel*param); % CPD energy function.
        dP=spdiags(sum(P,2),0,m,m); % precompute diag(P)
        
        % with ctrl_pts
        param =(basis'*dP*basis+lambda*sigma^2*kernel)\(basis'*(P*scene-dP*model0));
        
        
        % update model
        model = model0 + basis*param;
        EMtol = norm((E_old-E)/E_old);
        EMiter = EMiter + 1;
    end  % end of iteration/perT
    % Anneal
    sigma = sigma * anneal_rate;
    iter = iter + 1;
    hold on
    plot(model(:,1),model(:,2),'.')
    
    
    ntol = norm(model_old - model);
end % end of annealing.
toc
model = cpd_denormalize(model, centroid, scale);

end

%%
function [P, E] = cpd_P(x, y, sigma, outliers)

if nargin<3, error('cpd_P.m error! Not enough input parameters.'); end;
if ~exist('outliers','var') || isempty(outliers), outliers = 0; end;

k=-2*sigma^2;
[n, d]=size(x);[m, d]=size(y);

P=repmat(x,[1 1 m])-permute(repmat(y,[1 1 n]),[3 2 1]);

P=squeeze(sum(P.^2,2));
P=P/k;
P=exp(P);

% compute column sums -> s
if outliers
    Pn=outliers*(-k*pi)^(0.5*d)*ones(1,m);
    s=sum([P;Pn]);
else
    s=sum(P);
end

if nnz(s)==numel(s)
    E=-sum(log(s));  % log-likelihood
    P=P./repmat(s,n,1); % normalization such that each column sums to 1
    %    s=sum(P,2);
    %    P=P./repmat(s,1,m); % normalization such that each row sums to 1
else
    P=[];E=[];
end

end
%%
function G=cpd_G(x,y,beta)

if nargin<3, error('cpd_G.m error! Not enough input parameters.'); end;

k=-2*beta^2;
[n, d]=size(x); [m, d]=size(y);

G=repmat(x,[1 1 m])-permute(repmat(y,[1 1 n]),[3 2 1]);
G=squeeze(sum(G.^2,2));
G=G/k;
G=exp(G);

end

%%
function  [X, centroid, scale] = cpd_normalize(x)

[n, d]=size(x);
centroid=mean(x);
x=x-repmat(centroid,n,1);
scale=sqrt(sum(sum(x.^2,2))/n);
X=x/scale;
%X = x;
end
%%

function x =cpd_denormalize(X, centroid, scale)

[m, d]=size(X);
x=X*scale;         % scale
x=x+repmat(centroid,m,1); % move
end

%%  