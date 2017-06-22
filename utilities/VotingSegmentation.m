function [label_img,NumOfClusters] = VotingSegmentation(H,W,Nsegs,para,img_size,featPath)
%% This function apply the voting cut on the given bipartite graph
% Input:
% H - the incident matrix of the bipartite graph, a m x n binary matrix;
% W - the weights on the edges, a n x n diagonal matrix;
% Nsegs -  the number of clusters in the superpixel set, i.e. the dimension
% of the superpixel embedding space.
% img_size -  the image size
% Output:
% label_img - the hierarchical segmentation of the image, i.e. a cell, each element is a segmentation of the image

%% get superpixel associations
[H_r,~,ic] = unique(H,'rows','stable');
RegionSize = hist(ic,unique(ic))';
SpAssocSize = bsxfun(@times, H_r,RegionSize);
%NumOfRegions = length(RegionSize);

% IdxOfKeptRegion = RegionSize>para.minSize;
% KeptH_r = H_r(IdxOfKeptRegion,:);
% IdxOfRemovedRegion = IdxOfKeptRegion==0;
% RemovedH_r = H_r(IdxOfRemovedRegion,:);
% simMatrix = RemovedH_r*H_r';
% simMatrix = bsxfun(@times, simMatrix, IdxOfKeptRegion');
% [~,IdxOfMergedTo] =max(simMatrix,[],2); % The removedH_r will be merged into a keptH_r that is most similar to it.

%% embeding superpixel with spectral clustering (Normalized cut) 
B = sparse(SpAssocSize)*W;
[Nx,Ny] = size(B);
%%% build the superpixel graph
dx = sum(B,2);
Dx = sparse(1:Nx,1:Nx,1./dx);
% Dx = sparse(1:Nx,1:Nx,1./sqrt(dx));
% normB = Dx*B;
Wy = B'*Dx*B;
d = sum(Wy,2);
% D = sparse(1:Ny,1:Ny,1./d);
% normWy = D*Wy;
D = sparse(1:Ny,1:Ny,1./sqrt(d));
nWy = D*Wy*D;
nWy = (nWy+nWy')/2;


%%  voting cut

switch para.votingmode
    case 'VotingTransfer'
        % computer eigenvectors
        [evec,~] = eigs(full(nWy),Nsegs); % use eigs for large superpixel graphs  
        %Ncut_evec = D*evec(:,idx(1:Nseg));

        % compute the segmentation of superpixels
        % normalize each row to unit norm
        Norm_evec = bsxfun( @rdivide, evec, sqrt(sum(evec.*evec,2)) + 1e-10 );

        % k-means
        splabels = k_means(Norm_evec',Nsegs);
        % convert the labels into a binary matrix
        BinarySpLabels = zeros(length(splabels),Nsegs);
        for i = 1:Nsegs
            BinarySpLabels(:,i)=splabels==i;
        end 
        [~,InitialRegionLabels] = max(B*BinarySpLabels,[],2); % voting to get the superpixel associations' clustering label
        % transfer the segmentation information by voting 
        pixelLabels = InitialRegionLabels(ic); % convert the superpixel association label to pixels
        label_img = reshape(pixelLabels,img_size); % get the label mask 
        NumOfClusters = Nsegs;
    case 'HierarchicalTree'
        % merge the tiny superpixel association
        temp = unique(RegionSize);
        minSize = floor(quantile(temp,para.minSizeThreshold));
        tempRegionSize = RegionSize;
        tempThreshold = 1;
        while sum(tempRegionSize<=minSize)>8100 % for memory issue,lanch a quick merge to merge the tiny sp assoc into the nearest one. Practically, the memory threshold is set to be 8100. 
            fprintf('tempThreshold: %d\n',tempThreshold);
            if tempThreshold>1
                tempRegionSize = RegionSize;
                H_r = unique(H,'rows','stable'); % It is time consuming!!!
            end
            IdxOfKeptRegion = tempRegionSize>tempThreshold;% set the minSize to 10.
            %KeptH_r = H_r(IdxOfKeptRegion,:);
            IdxOfRemovedRegion = IdxOfKeptRegion==0;
            RemovedH_r = H_r(IdxOfRemovedRegion,:);
            simMatrix = RemovedH_r*H_r';
            simMatrix = bsxfun(@times, simMatrix, IdxOfKeptRegion');
            [~,IdxOfMergedTo] =max(simMatrix,[],2); % The removedH_r will be merged into a keptH_r that is most similar to it.
            H_r(IdxOfRemovedRegion,:)=H_r(IdxOfMergedTo,:);
            H_rebuilt = H_r(ic,:);
            [H_r,~,ic] = unique(H_rebuilt,'rows','stable'); 
            tempRegionSize = hist(ic,unique(ic))';
            tempThreshold = tempThreshold+1;
        end
        RegionSize = tempRegionSize;
%         while sum(RegionSize<=minSize)>8100 % for memory issue,lanch a quick merge to merge the tiny sp assoc into the nearest one. 
%         
%         end
        %fprintf('the minSize is %d\n',minSize);
        [NewH,NewSize,Parent,Kept] = mergingSpAssoc(H_r,RegionSize,minSize);
        % computer eigenvectors
        [evec,~] = eigs(full(nWy),max(para.rangeOfK)); % use eigs for large superpixel graphs  
        % Voting for segmentation 
        %Tree = HierarchicalSegTreeVer2(evec,para.rangeOfK,KeptH_r);
        Tree = HierarchicalSegTreeVer2(evec,para.rangeOfK,NewH,NewSize);
        %[label_img,NumOfClusters]=getLabel(Tree,ic,IdxOfKeptRegion,IdxOfRemovedRegion,IdxOfMergedTo,img_size);
        [label_img,NumOfClusters] = getLabel(Tree,ic,Parent,Kept,img_size,para.lifetimelevel);
%         label_img = cell(1,length(Tree));
%         for i = 1: length(Tree)
%             label = Tree{i}(ic);
%             label_img{i}=reshape(label,img_size);
%         end
        % label_img = reshape(pixLabels(1:img_size(1)*img_size(2)),img_size);
end
end



