function [conF] = Tuning(conF,vecPosCon)
% Ajuste fino de trayactori
for jj = 1:numel(conF)
    if jj == numel(conF)
%         Co1= Translation_matrix(vecPosCon{end,2}',[ 0 0 0])';
%         Co2 = Translation_matrix(conF{jj}{end}',[ 0 0 0])';
        Co1 = vecPosCon{end,2}';
        Co2 = conF{jj}{end};
        [~,~,T,~]=icp(Co1([3 1],:),Co2([3 1],:)); % vector de ajuste
        t = linspace(0,1,numel(conF{jj,1}));
        for ii = 1:numel(conF{jj,1})
            conF{jj,1}{ii,1}=Translation_matrix((conF{jj,1}{ii,1}+[T(1) 0 T(2)]*t(ii))',mean(conF{jj,1}{ii,1}))'; % ajustar solo en X and Y
        end
    else
        Co1= Translation_matrix(conF{jj+1}{1}',[ 0 0 0])';
        Co2 = Translation_matrix(conF{jj}{end}',[ 0 0 0])';
        
        [~,~,T,~]=icp(Co1([3 1],:),Co2([3 1],:)); % vector de ajuste
        t = linspace(0,1,numel(conF{jj,1}));
        for ii = 1:numel(conF{jj,1})
            
            conF{jj,1}{ii,1}= Translation_matrix((conF{jj,1}{ii,1}+[T(1) 0 T(2)]*t(ii))',mean(conF{jj,1}{ii,1}))'; % ajustar solo en X and Y
        end
    end
end