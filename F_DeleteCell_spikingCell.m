% clear all
% close all

% k=1;
% run('C:\UoM\Adrien Data\AllData\Load_Behavior.m'); 
% ITI = load(['NewITI_NoReward_',num2str(k),'.txt']);
% TrData = cell(size(Behav{k},1)-1,1);
% for tr = 1:size(Behav{k},1)-1
%     indITI = find((ITI(:,2)>=Behav{k}(tr,2)) & (ITI(:,2)<Behav{k}(tr+1,1)));
%     TrData{tr} = ITI(indITI,:);
% end
% run('Load_CellType.m'); 
% celltipo = CellType{k};
% clear CellType
% CellType = celltipo;
% thresh = 0;
% SpkData = TrData;
% nCell = length(CellType);


function [SpkData, num, nallIDs] = F_DeleteCell_spikingCell(SpkData,nCell,thresh)
% 
% This function delete the cells that are silent for more than a specific 
% percentage (thresh) of trials along the session 
% INPUT
% SpkData :   are the spike trains of each cell for each trial. This is
%             a cell array. For each trial a cell contains the spike 
%             trains of the neurons   
% nCell   :   is a scalar and contains the total number of cell recorded in
%             the session
% thresh  :   is the threshold between 0 and 1 over which the silent cell 
%             will be deleted. 0 if a cell is silent in at least one trial 
%             is eliminated; 1 no cell are eliminated. 
% OUTPUT
% SpkData :   new cell array containing only the spike trains that satisfy
%             the treshold condition
% num     :   is a vector containing the index of the cell excluded


    emptyspk = [];
%     nall = zeros(size(SpkData,1),1);
    for j = 1:size(SpkData,1)
        SpikeData = SpkData{j};

        IDspk = unique(SpikeData(:,1));
%         nall(j) = numel(IDspk);
        nallIDs = numel(IDspk);
        if nCell > nallIDs
           ISM = ismember(1:nCell,IDspk);
           emptyspk = [emptyspk find(ISM==0)];
        end
        
        for i = 1:nCell
            if sum(SpikeData(:,1)==i)==1
                emptyspk = [emptyspk i];
            end
        end
                
    end
%     nallIDs = min(nall);
    
    IDemptyspk = unique(emptyspk);
    num = [];
    for h = 1:length(IDemptyspk)
        if length(find(emptyspk==IDemptyspk(h)))/size(SpkData,1) > thresh
            num(end+1) = IDemptyspk(h);
%             indspk = IDemptyspk(h);
            for j = 1:size(SpkData,1)
                clear SpkData{j}
                SpkData{j}(SpkData{j}(:,1)==num(end),:) = [];
            end
        end
    end
    
end
    
    
    
 