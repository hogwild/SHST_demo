function [label_img, clusterNumber] = getLabel(SegTree,ic,parentId,theKeptId,img_size,lifetimelevel)%IdxOfKeptRegion,IdxOfRemovedRegion,IdxOfMergedTo,img_size,threshold) 
%% This function is label the pixels based on the SegTree
% Input:
% SegTree - A hierarchical segmentation tree (or forest).Note:it is a maximum spanning
% tree
% ic - the (pixel) - (superpixel association) map
% img_size - size of the image
% threshold - the threshod for partition the tree;default is the cluster number of the highest lifetime.
% Output:
% label - the label of the pixels, an matrix 
% NumberOfClusters - the cluster number
%% get labels for superpixel association
import java.util.LinkedList;
[~,~,edgeWeights]=find(SegTree);
edgeWeights = sort(edgeWeights);
lifetime = diff(edgeWeights);
sortedlifetime = sort(unique(lifetime),'descend');
if nargin<6
    selectlifetime = sortedlifetime(1);% i.e. the highest lifetime
    %[~,idx] = max(lifetime);
    threshold = edgeWeights(lifetime==selectlifetime);
    threshold = max(threshold);% to avoid there are more than one threshold;
else
    selectlifetime = sortedlifetime(lifetimelevel);
    threshold = edgeWeights(lifetime==selectlifetime);
    threshold = max(threshold);% to avoid there are more than one threshold;
end

% pick out the root of subtree
[rows,cols] = find((SegTree>0)&(SegTree<=threshold)); % note that we use a max spanning tree, so should be '<=' and 0 should not be considered as a edge value.
% remove the edges between the subtrees
SegTree(rows,cols)=0; 
subTreeRoot = getRoot(SegTree);
%subTreeRoot = [rows;cols];
%subTreeRoot = unique(subTreeRoot);
% preorder tranversal (DLR?) of each subtree;
label = zeros(size(SegTree,1),1);
queueSubTreeRoot = LinkedList();
queue = LinkedList();
clusterNumber = 0;
for i = 1:length(subTreeRoot)
    queueSubTreeRoot.offer(subTreeRoot(i));
end
%ShouldRemove = -1;
%Parent = getParent(SegTree);
while ~isempty(queueSubTreeRoot.peek())
    clusterNumber = clusterNumber + 1;
    currentRoot = queueSubTreeRoot.poll();
    queue.offer(currentRoot);
    while ~isempty(queue.peek())
        current = queue.poll();
        label(current) = clusterNumber;
        children = getChildren(SegTree,current);
        for child = 1:length(children)
            queue.offer(children(child));
        end
    end
    
end
label_spAssoc = zeros(length(theKeptId),1);
label_spAssoc(theKeptId)=label;
label_spAssoc =label_spAssoc(parentId);
label_img = label_spAssoc(ic);
label_img = reshape(label_img,img_size);
% label_spAssoc = zeros(size(IdxOfKeptRegion,1),1);
% label_spAssoc(IdxOfKeptRegion) = label;
% label_spAssoc(IdxOfRemovedRegion) = label_spAssoc(IdxOfMergedTo);
% label_img = label_spAssoc(ic);
% label_img = reshape(label_img,img_size);
%clear queueSubTreeRoot;
%clear queue;

end

%% sub function: get the children of the given vertex
function children = getChildren(Tree,idx)
children = find(Tree(:,idx)~=0);
end

%% sub function: get the root of the given tree matrix(note there may be more than one root,i.e a few of sub trees)
function [root] = getRoot(Tree)
root = find(sum(Tree,2)==0);
end
