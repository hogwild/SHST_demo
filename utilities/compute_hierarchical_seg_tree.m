function  [SegTree,Parent,Kept] = compute_hierarchical_seg_tree(eigvec,RangeOfK,H_r,RegionSize,minSize)
%% this function is compute the hierarchical tree for the superpixel associations
% Input:
% eigvec - eigenvectors of the superpixel
% RangeOfK - the range of cluster number k
% H_r = the superpixel association
% Output:
% SegTree - The hierarchical segmentation tree
% Lifetime - the longest lifetime, for cluster number decision


%% parameters
NumOfK = length(RangeOfK);
NumOfSpAssoc = size(H_r,1);
SpAssocHistOfSp = bsxfun(@times,H_r,RegionSize) ;
%NumOfVoters = sum(H_r,2);


%CoAssocMatrix = zeros(NumOfSpAssoc);
SpAssocHistAll = [];

%% compute the co-association matrix
for i = 1: NumOfK;
    k = RangeOfK(i);
    Norm_evec = bsxfun( @rdivide, eigvec(:,1:k), sqrt(sum(eigvec(:,1:k).*eigvec(:,1:k),2)) + 1e-10 );
    splabels = k_means(Norm_evec',k);
    SpAssocClusterIdx = bsxfun(@times, H_r, splabels');
    SpAssocHist = zeros(NumOfSpAssoc,k);
    for j = 1:k
        votes = SpAssocClusterIdx == j;
        bin = sum(SpAssocHistOfSp.*votes,2);%./NewSize; %normalized bin value
        SpAssocHist(:,j)=bin;
    end
    SpAssocHistAll = [SpAssocHistAll,SpAssocHist];
end
% merge the tiny sp assoc
[MergedSpAssocHistAll,NewSize,Parent,Kept] = merging_SpAssociation(SpAssocHistAll,RegionSize,minSize);
% compute the co-association matrix
CoAssocMatrix = SpAssoc_SimMeasure(MergedSpAssocHistAll,NewSize);

%% compute the Minimum Span Tree
% check the connection of the graph
[NumOfTrees,root] = Check_Graph_Connection(CoAssocMatrix);
if NumOfTrees == 1
    [SegTree,~] = mst(-CoAssocMatrix);
else
    SegTree = zeros(size(CoAssocMatrix));
    for i = 1:length(root)
        [SegTree_i,~] = mst(-CoAssocMatrix,root(i)); % in fact we need a maximum spanning tree, because we use minimum span tree, we convert it to a minus one.
        SegTree = SegTree + SegTree_i;
    end
end
SegTree = -SegTree; 

%% get labels for superpixel association
% [~,~,edgeWeights]=find(SegTree);
% edgeWeights = sort(edgeWeights);
% lifetime = diff(edgeWeights);
% highestlifetime = max(lifetime);
% % pick out the root of subtree
% [rows,cols] = find(SegTree>highestlifetime);
% % remove the edges between the subtrees
% SegTree(rows,cols)=0; 
% subTreeRoot = getRoot(SegTree);
% %subTreeRoot = [rows;cols];
% %subTreeRoot = unique(subTreeRoot);
% % preorder tranversal (DLR?) of each subtree;
% label = zeros(size(SegTree,1),1);
% queueSubTreeRoot = LinkedList();
% queue = LinkedList();
% clusterNumber = 0;
% for i = 1:length(subTreeRoot)
%     queueSubTreeRoot.offer(subTreeRoot(i));
% end
% %ShouldRemove = -1;
% %Parent = getParent(SegTree);
% while ~isempty(queueSubTreeRoot.peek())
%     clusterNumber = clusterNumber + 1;
%     currentRoot = queueSubTreeRoot.poll();
%     queue.offer(currentRoot);
%     while ~isempty(queue.peek())
%         current = queue.poll();
%         label(current) = clusterNumber;
%         children = getChildren(SegTree,current);
%         %parent = Parent(current);
%         %connectedVertex = [children;parent];
% %         for i = 1:length(subTreeRoot) % check if the vertex in children exists in subTree root.If it exists, then remove it from the subTree root queue.
% %             if sum(connectedVertex == subTreeRoot(i))>0
% %                queueSubTreeRoot.remove(i-clusterNumber);
% %                ShouldRemove = i; 
% %             end
% %         end
% %         if ShouldRemove~=-1
% %             subTreeRoot(ShouldRemove)=[];
% %             ShouldRemove = -1;
% %         end
%         for child = 1:length(children)
%            % if label(connectedVertex(i))==0 % to check if the vetex has been labeled.If it has been labeled, ignore it.
%                 queue.offer(children(child));
%             %end
%         end
%     end
%     
% end
% %clear queueSubTreeRoot;
% %clear queue;
% end
% 
% %% sub function: get the children of the given vertex
% function children = getChildren(Tree,idx)
% children = find(Tree(:,idx)~=0);
% end
% 
% %% sub function: get the root of the given tree matrix(note there may be more than one root,i.e a few of sub trees)
% function [root] = getRoot(Tree)
% root = find(sum(Tree,2)==0);
% end







