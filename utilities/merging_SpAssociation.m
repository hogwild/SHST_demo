function [SpAssocHistAll, NewSize, ParentId,theKeptId] = merging_SpAssociation(SpAssocHistAll,SpAssoc_Size,minSize)
%% this function merge the superpixel association whose size is small than minSize into the most similar one
% Input:
% H_r: the superpixel association
% SpAssoc_Size: the size of every superpixel association
% minSize: the size threshold
% Output:
% New_H: the merged superpixel association list
% IdMaps: a vector of the id of the superpixel association in the old H_r.
% so, label_of_H_r_rows = label_of_New_H_rows(IdMaps).
%% preparing
import java.util.LinkedList;
NumOfSpAssoc = size(SpAssocHistAll,1);
removed_Idx = SpAssoc_Size<=minSize;
RemovedH = SpAssocHistAll(removed_Idx,:);
RemovedIdx = find(removed_Idx);
queueRemovedH = LinkedList();
queueRemovedIdx = LinkedList();
ParentId = 1:NumOfSpAssoc;
% if sum(SpAssoc_Size<=minSize)>8100 %% for memory issue,lanch a quick merge to merge the tiny sp assoc into the nearest one. 
%     removed_Idx = SpAssoc_Size<10; %% set the minSize to 10
%     
% 
% end

for i = 1:size(RemovedH,1);
    queueRemovedH.offer(RemovedH(i,:));
    queueRemovedIdx.offer(RemovedIdx(i));
end
%kept_Idx = removed_Idx==0;
%NewH = H_r(kept_Idx);

while ~isempty(queueRemovedH.peek())
    currentH = queueRemovedH.poll(); % here, the currentH is a column vector
    currentIdx = queueRemovedIdx.poll();
    if isequal(currentH',SpAssocHistAll(currentIdx,:))% to check if currentH has absorbed other spAssoc.
        SpAssocHistAll(currentIdx,:)=0; % set the currentH = 0, i.e. remove it but unchange the size of NewH.
        sim = simMeasure(SpAssocHistAll,currentH,SpAssoc_Size,SpAssoc_Size(currentIdx)); % get the similarity to other H.
        [~,mergeToIdx] = max(sim);
        ParentId(ParentId==ParentId(currentIdx)) = mergeToIdx; %change the parent id
        SpAssocHistAll(mergeToIdx,:) = currentH' + SpAssocHistAll(mergeToIdx,:);
        SpAssoc_Size(mergeToIdx) = SpAssoc_Size(mergeToIdx) + SpAssoc_Size(currentIdx);
        %SpAssoc_Size(currentIdx) = 0;% set the current SpAssoc size to 0.
        if SpAssoc_Size(mergeToIdx)<=minSize
% unable to locate the tiny sp assoc,so , use an if statement at the begining of while instead.
%             queueRemovedH.remove(?); 
%             queueRemoveIdx.remove(?);
            queueRemovedH.offer(SpAssocHistAll(mergeToIdx,:));
            queueRemovedIdx.offer(mergeToIdx);
        end
    end
end
theKeptId = sum(SpAssocHistAll,2)>0; 
SpAssocHistAll = SpAssocHistAll(theKeptId,:);
NewSize = SpAssoc_Size(theKeptId);
end

%% sub function for similarity measure
function sim = simMeasure(HistArray,currentHist,RegionSize,CurrentRegionSize)
HistInnerProduct = HistArray*currentHist;
Size = RegionSize*CurrentRegionSize;
sim = HistInnerProduct./Size;
end

