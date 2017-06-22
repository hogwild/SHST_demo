function dist = chiSq_dist(histMat)
%% Compute the chi-square distance between the rows of a given histgram matrix
% Input:
% histMat - a (m x n) histgram matrix, each row is a histgram with n bins
% Output:
% dist = a (m x m) matrix, the i-th row is the distance of node i to the other nodes.
%%
[m,~]=size(histMat);
% a fast speed but high memory cost version
% if issparse(histMat)
%     histMat=full(histMat);
% end

% a low speed but memory saving version
% if issparse(histMat)~=1
%     histMat=sparse(histMat);
% end
total = sum(histMat,2);
total(total==0)=1;
histMat=bsxfun(@rdivide,histMat,total);
dist = zeros(m,m);
for i = 1:m
    temp1 = bsxfun(@minus, histMat,histMat(i,:));
    temp2 = bsxfun(@plus, histMat,histMat(i,:));
    temp2(temp2==0)=1; % set those bin = 0 in both two histogram to 1. In this case temp1 also equal to 0, so, the result temp1.^0.5./temp2==0
    temp3 = (temp1.^2)./temp2;
    dist(i,:) = sum(temp3,2);
end
%dist = dist.*0.5;

