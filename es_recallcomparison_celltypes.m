close all
clearvars

ndata = 50;
learnSes = [1 3 7 14 19 28 32 41 44 46];
% InitSes = [1 13 27 38];

% type = 'Trial';

% Load recall data for all cells and only pyramidal
deltarec_all = csvread('Trial_ErrorCorrect_pvalue_difference_allcells.csv');
deltarec_pyr = csvread('Trial_ErrorCorrect_pvalue_difference_onlypyramidal.csv');

% visualize the recall results by comparing the two dataset
jit = randn(ndata,1)*0.05;
figure()
h1 = plot([ones(ndata,1)+jit ones(ndata,1)*2+jit]',[deltarec_all(:,3) deltarec_pyr(:,3)]',...
    'o-','Color',[.7 .7 .7]); hold on
h2 = plot([ones(length(learnSes),1)+jit(learnSes) ones(length(learnSes),1)*2+jit(learnSes)]',[deltarec_all(learnSes,3) deltarec_pyr(learnSes,3)]',...
    'ro-'); hold on
plot([.5 2.5],[0 0],'k--')
legend([h1(1), h2(1)],'Other','Learning')
set(gca,'XTick',1:2)
set(gca,'XTickLabel',{'All cells','Pyramidal'})
ylabel('Delta Recall')


% insert stats
[h,p]= ttest2(deltarec_all(:,3), deltarec_pyr(:,3))
[h,p]= kstest2(deltarec_all(:,3), deltarec_pyr(:,3))
[p,h]= ranksum(deltarec_all(:,3), deltarec_pyr(:,3))

[h_alllearn,p_alllearn]= kstest(deltarec_all(learnSes,3))
[h_pyrlearn,p_pyrlearn]= kstest(deltarec_pyr(learnSes,3))
[h_all,p_all]= kstest(deltarec_all(:,3))
[h_pyr,p_pyr]= kstest(deltarec_pyr(:,3))