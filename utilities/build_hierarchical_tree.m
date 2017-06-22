function [Tree,Parent,Kept,ic] = build_hierarchical_tree(H,B,para)
%% This function apply the voting cut on the given bipartite graph
% Input:
% H - the incident matrix of the bipartite graph, a m x n binary matrix;
% B - the weighted adjacency matrix, for spectral clustering on superpixels
% para - the paremeter set;
% Output:
% Tree - a hierarchical tree;
% Parent - the parent of the removed nodes index in H_r
% Kept - the kept node index in H_r
% ic - the index from H_r to H
% Note: in this version, the eigenvectors of Ly is symetricly normalized, we can also convert them into the random walk normalized version. 
%% embeding superpixel with spectral clustering (Normalized cut) 
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
%% get superpixel associations
[H_r,~,ic] = unique(H,'rows','stable');
RegionSize = hist(ic,unique(ic))';
%% Hierarchical tree
% merge the tiny superpixel association
temp = unique(RegionSize);
minSize = floor(quantile(temp,para.minSizeThreshold));

% computer eigenvectors
[evec,~] = eigs(full(nWy),max(para.rangeOfK)); % use eigs for large superpixel graphs  
% evec = D*evec; % Note: by this step, the eigenvector will turn into a
% random walk normalized version.
% cut form, need to try in future (just make evec = D*evec). I guess the results won't change to much
% Voting for segmentation 
%Tree = HierarchicalSegTreeVer2(evec,para.rangeOfK,H_r,RegionSize);
[Tree,Parent,Kept] = compute_hierarchical_seg_tree(evec,para.rangeOfK,H_r,RegionSize,minSize);
        
end
