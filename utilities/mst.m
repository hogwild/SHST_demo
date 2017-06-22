function [Tree,parents] = mst(G,root)
% The minimum spanning tree algorithm (using the Prim's algorithm).
% Yangqing Jia, jiayq@eecs.berkeley.edu, 2010

if nargin==1
    % if root is not defined, simply let the first vertex to be the root -
    % it does not make any difference since picking any vertex in the mst
    % yields to the same set of edges.
    root = 1;
end

inset = false(1,size(G,1));
outset = true(1,size(G,1));
inset(root) = true;
outset(root) = false;

Tree = zeros(size(G));
parents = zeros(1,size(G,1));
for iter = 2:size(G,1)
    [ignore,idxi,idxj] = findmin(G,inset,outset);
    Tree(idxj,idxi) = G(idxj,idxi);
    parents(idxj) = idxi;
    inset(idxj) = true;
    outset(idxj) = false;
end

end
%% subfunction: findmin
function [minval,idxi,idxj] = findmin(G,rows,columns)
[i,j,val] = find(G(rows,columns));
[minval,idx] = min(val);
idxrows = find(rows);
idxcols = find(columns);
idxi = idxrows(i(idx));
idxj = idxcols(j(idx));
end