%% Load data
tmp=load('p1_RD_twin_database.mat')
grains_p1_RD=tmp.grains;
stats_p1_RD=tmp.stats;
twin_p1_RD=tmp.twin; 
mergedGrains_p1_RD=tmp.mergedGrains; 

tmp=load('p1_TD_twin_database.mat')
grains_p1_TD=tmp.grains;
stats_p1_TD=tmp.stats;
twin_p1_TD=tmp.twin; 
mergedGrains_p1_TD=tmp.mergedGrains; 

tmp=load('p1_ND_twin_database.mat')
grains_p1_ND=tmp.grains;
stats_p1_ND=tmp.stats;
twin_p1_ND=tmp.twin; 
mergedGrains_p1_ND=tmp.mergedGrains; 
%%
% hist(grains_p1_TD.prop.schmidActiveRank
% grains_p1_TD.prop.schmid(:,grains_p1_TD.prop.schmid)
% bpcombined = [bpcombine1(:), bpcombine2(:), bpcombine3(:)];
% hb = bar(xdata, bpcombined, 'grouped')
% 
% data = [randi([1 3], 5, 1)  randi([4 6], 5, 1)  randi([4 6], 5, 1)];
% figure(1)
% hb = bar(data)
% set(hb(1), 'FaceColor','r')
% set(hb(2), 'FaceColor','b')
% set(hb(3), 'FaceColor','g')