clear; close all; clc
load mm.mat

pca(conFT{end})
pca(vecPosCon{2,2}) 


 [TR,TT,TT2,data] = icp(conFT{end}(:,[1 3]),vecPosCon{2,2}(:,[1 3]),100,5,4);
figure(33)
hold on
scatter(conFT{end}(:,1),conFT{end}(:,3),'+')
scatter(vecPosCon{2,2}(:,1),vecPosCon{2,2}(:,3),'*')
scatter(data(1,:),data(2,:),'O')