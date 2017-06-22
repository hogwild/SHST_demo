function A = prune_knn(W,nb)

% W - table of similarity in size of n x n (may not symetric), the i_th row represent the
% similarity of the i_th superpixel to the others. 
% nb - the number of nearest neighbors 
% A - is a sparse matrix. represent the weights on the edges.

n = size(W,1);

assert(nb<n,'error:Too much neighbours!');

[~,idx] = sort(W,2,'descend');

idxC = idx(:,1:nb);
idxR = repmat((1:n)',1,nb);%ones(nb,1)*[1:n];
A = sparse(idxR(:),idxC(:),1,n,n); % for the i_th row(represent the i_th node), locates the nb_est similar nodes of the i_th node.
A = double(A|A'); %make the out symetric.
A = A.*W;

