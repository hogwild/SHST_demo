function label_img = Vcut(B,Nseg,img_size)

% B - |X|-by-|Y|, cross-affinity-matrix
% note that |X| = |Y| + |I|

[Nx,Ny] = size(B);
if Ny < Nseg
    error('Need more superpixels!');
end

%%% build the superpixel graph
dx = sum(B,2);
Dx = sparse(1:Nx,1:Nx,1./dx);

%NormalizedB=sparse(1:Nx,1:Nx,1./sqrt(dx))*B;
%D_e = sparse(inv(diag(sqrt(sum(B>0)))));
%NormalizedB=NormalizedB*D_e;
Wy = B'*Dx*B;
%Wy = NormalizedB'*NormalizedB;
%%% compute Ncut eigenvectors
% normalized affinity matrix
d = sum(Wy,2);
D = sparse(1:Ny,1:Ny,1./d);
normWy = D*Wy;
D = sparse(1:Ny,1:Ny,1./sqrt(d));
nWy = D*Wy*D;
nWy = (nWy+nWy')/2;

% computer eigenvectors
% [evec,~] = eigs(full(nWy),2);
[evec,~] = eigs(full(nWy),Nseg); % use eigs for large superpixel graphs  
%Ncut_evec = D*evec(:,idx(1:Nseg));

% compute the segmentation of superpixels



%%% compute the Ncut eigenvectors on the entire bipartite graph (transfer!)
%evec = Dx * B * Ncut_evec;

%%% k-means clustering
% extract spectral representations for pixels
%evec = evec(1:prod(img_size),:);

% normalize each row to unit norm
evec = bsxfun( @rdivide, evec, sqrt(sum(evec.*evec,2)) + 1e-10 );

% k-means
splabels = k_means(evec',Nseg);
BinarySpLabels = zeros(length(splabels),Nseg);
for i = 1:Nseg
    BinarySpLabels(:,i)=splabels==i;
end 
% Voting for segmentation 
[~,pixLabels] = max((B>0)*normWy*BinarySpLabels,[],2); % the vote transfer
%[~,pixLabels] = max((B>0)*nWy*evec,[],2);
label_img = reshape(pixLabels(1:img_size(1)*img_size(2)),img_size);

