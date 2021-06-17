clear all 
% close all

run('Load_CellType.m') 
run('Load_Behavior.m') 

Qt = 0.01; % (seconds) time-resolution for convolution window : here 10 ms
% cluster analysis function arguments
binlessopts.Dmeth = 'corr';
binlessopts.BLmeth = 'Gaussian';
binlessopts.modopts = {{'sqEuclidean'},100};  % use 100 repetitions of k-means with Euclidean distance as basis for consensus
binlessopts.BLpars = 0.1; % width of convolution window, in seconds (here SD of Gaussian)
binlessopts.blnS = 0; % if 1 does not do the clustering, but compute only computes the similarity matrix
startts = 0;    % proportion of time from start of recording to begin analysis [0=first spike]
endts = 1;      % proportion of time from end of recording to end analysis [1 = final spike]

ndata = 50;
thresh = 0;
deconstructThresh = .7;
start = 1;
learnSes = [1 3 7 14 19 28 32 41 44 46];
InitSes = [1 13 27 38];

learning = zeros(ndata,1);

num = cell(ndata,1);
nallIDs = cell(ndata,1);
Rmatrix = cell(ndata,1);
before = cell(ndata,1);
after = cell(ndata,1);
newBehav = cell(ndata,1);
Rmatrice = cell(ndata,1);
type = 'Trial';

for k = 1:ndata

    switch type
        case 'ITI'
            load(['Processed Data\ITI_SpikeTrainsXTrial_',num2str(k)])
            TrData = ITIData;
        case 'Trial'
            load(['Processed Data\SpikeTrainsXTrial_',num2str(k)])
    end
    
    Cxy = cell(size(TrData,1),1);
    Sxy = cell(size(TrData,1),1);
    
    [TrData, num{k}, nallIDs{k}] = F_DeleteCell_spikingCell(TrData,length(CellType{k}),thresh);
    
%      %% select the interneurons (1) or pyramidal (2) neurons active in every trials or ITIs
%     pyr = find(CellType{k}==2);
%     CellIndex = setdiff(pyr,num{k});
    
    %% select all the neurons active in every trials or ITIs 
    CellIndex = setdiff(1:length(CellType{k}),num{k});
    
    IDspk = CellIndex;    
    N = length(IDspk);
    
    SpikeTrain = cell(size(TrData,1),1);
    NewTrData = cell(size(TrData,1),1);
    for j = 1:size(TrData,1)
        
        % Z-score of data: In this case I have used the mean and std over
        % all spike trains.
        TrData{j}(:,2) = (TrData{j}(:,2)-mean(TrData{j}(:,2)))./sqrt(var(TrData{j}(:,2)));
        
        T_start_recording = min(TrData{j}(:,2));
        T_end_recording = max(TrData{j}(:,2));
        T_period = T_end_recording - T_start_recording;
        
        % fix start and end times of data-set to use
        T = [T_start_recording + startts*T_period T_start_recording + endts*T_period];

        %%%%%% now not-rectify similarity matrix %%%%%%%%%%%%%%%%%%%
        if ~exist('TrData','var'), error(['Data-file ' TrData{j} ' does not contain spks variable']); end  
        
        if ~isempty(IDspk)
            [Cxy{j},spkfcn,shiftbase] = F_SimilarityMatrix(TrData{j},IDspk,T,Qt,binlessopts);
            % rectify the similarity matrix
            Sxy{j}{1} = Cxy{j}{1};
            Sxy{j}{1}(Sxy{j}{1} < 0) = 0; 
            Sxy{j}{1}(eye(N)==1) = 0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    newIDspk{k} = IDspk;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % ***************************************
    switch type
        case 'ITI'
            erID = find(Behav{k}(1:end-1,4)==0);
            corID = find(Behav{k}(1:end-1,4)==1);
        case 'Trial'
            erID = find(Behav{k}(:,4)==0);
            corID = find(Behav{k}(:,4)==1);
    end
    learning(k) = length(erID)+1;
    newBehav{k} = Behav{k}([erID; corID],:);
    % ***************************************

    
    Rmatrix{k} = zeros(size(TrData,1),size(TrData,1));
    if ~isempty(Sxy{j})
        for j = 1:size(Rmatrix{k},1)   
            for l = 1:size(Rmatrix{k},1)  
                Rmatrix{k}(j,l) = corr2(Sxy{j}{1},Sxy{l}{1});
            end
        end
        Rmatrix{k}(eye(size(TrData,1))==1) = 0;
    
    end
    clear Idx
    % ***************************************
    % Reshape Rmatrix according to error and correct trial
    Rmatrix{k} = Rmatrix{k}(:,[erID; corID]);
    Rmatrix{k} = Rmatrix{k}([erID; corID],:);    
    % ***************************************    
    group = [];
    group(1:(length(erID)-start)^2) = 1; 
    before{k} = reshape(Rmatrix{k}(start:length(erID),start:length(erID)),(length(erID))^2,1);
    after{k} = reshape(Rmatrix{k}(length(erID)+1:end,length(erID)+1:end),(size(Rmatrix{k},1)-length(erID))^2,1);
    group((length(erID)-start)^2+1:(length(erID)-start)^2+(size(Rmatrix{k},1)-length(erID)+1)^2) = 2;

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Some STATS
    [h(k), p(k)] = kstest2(before{k}, after{k});
    [pv(k), hi(k)] = ranksum(before{k}, after{k});
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

% save_fname = ['Processed Data\' type '_RetainedCell'];
% save(save_fname,'newIDspk')
    
% fname = 'Processed Data\Rmatrix_ErrorCorrectCorrelation_ITI_47Ses.mat';
% save(fname,'before','after')
%  
% fname = ['Processed Data\' type '_Rmatrix_correlation_withnegativevalues'];
% save(fname,'Rmatrice')

% matlabpool close

Prima = zeros(ndata,1);
Dopo = zeros(ndata,1);
ErrCor = zeros(ndata,3);

figure()
for k = 1:ndata
    clear erID
    erID = find(Behav{k}(1:end-1,4)==0);
    Prima(k) = nanmean(before{k});
    Dopo(k) = nanmean(after{k});
    ErrCor(k,:) = [k p(k) Dopo(k)-Prima(k)];
    
    subplot(1,4,sum(k>=InitSes))        
    if p(k)<0.05 && ismember(k,learnSes)
        plot([1 2],[nanmean(before{k}) nanmean(after{k})],'ro-','LineWidth',2);hold on
    elseif p(k)<0.05 && ~ismember(k,learnSes)
            plot([1 2],[nanmean(before{k}) nanmean(after{k})],'o-','Color',[0.7 .7 .7],'LineWidth',2);hold on
    elseif ismember(k,learnSes)
           plot([1 2],[nanmean(before{k}) nanmean(after{k})],'ro--','LineWidth',2);hold on
    else
        plot([1 2],[nanmean(before{k}) nanmean(after{k})],'o--','Color',[0.7 .7 .7],'LineWidth',2);hold on
    end
    
%     if ErrCor(k,3)<=0
%         plot([1 2],[nanmean(before{k}) nanmean(after{k})],'ro--','LineWidth',2);hold on
%     end
    xlim([0.5 2.5])
    set(gca,'XTick',1:2)
    set(gca,'XTickLabel',{'Error','Correct'})
    xtickangle(45)
    title(['Rat ',num2str(sum(k>=InitSes))])
    
end
subplot(1,4,1)  
ylabel('Average CC')

jit = randn(ndata,1)*0.05;
figure()
subplot(1,2,1)
h1 = plot([Prima Dopo]','o-','Color',[.7 .7 .7]); hold on
h2 = plot([Prima(learnSes) Dopo(learnSes)]','ro-'); hold on
legend([h1(1), h2(1)],'Other','Learning')
set(gca,'XTick',1:2)
set(gca,'XTickLabel',{'Error','Correct'})
ylabel('Recall')
xlim([.5 2.5])
title('a)                                                                        ')

subplot(1,2,2)
plot(ones(ndata,1)+jit,Dopo-Prima,'o','Color',[.7 .7 .7]); hold on
plot(ones(length(learnSes),1)+jit(learnSes),Dopo(learnSes)-Prima(learnSes),'ro'); hold on
plot([.5 1.5],[0 0],'k--')
ylabel('Delta Recall')
ylim([-.2 .4])
title('b)                                                                        ')


% csvwrite('ITI_ErrorCorrect_pvalue_difference_47Ses.csv',ErrCor)
csvwrite('Trial_ErrorCorrect_pvalue_difference_allcells.csv',ErrCor)
% csvwrite('Trial_ErrorCorrect_pvalue_difference_onlypyramidal.csv',ErrCor)
