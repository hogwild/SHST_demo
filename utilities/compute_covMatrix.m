function c = compute_covMatrix(X,seg_k)
% X is a set of random variable, whose size is m*n*k. k is the number of
% the random variables, and m*n is the observations of Xk.
% seg_k is the group index of those Xi.
% the output c is the covariance matrix of those group Xi.
[m,n,k]=size(X);
X_reshape = reshape(X,m*n,k);
Nsp=length(seg_k);
c=cell(Nsp,1);
    for i = 1:Nsp  
        Indx=seg_k{i};
        c{i}=cov(X_reshape(Indx,:));
    end
end