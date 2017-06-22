function distMat = compute_ChiSqDist(texton,seg_k)
%% this function compute the distance between superpixels base on texton feature
% Input:
% texton - texton feature of every pixel, it is an array in size of
% m*n( m*n is the size of the image)
% seg - the superpixel labels of the image.
% Output:
% distMat -  a symmetric matrix of distance between superpixels.
%% get the histogram of texton of each superpixel
Nsp=length(seg_k);
NumOfBins = max(texton(:))+1; % note the texton is start from 0, so, the bin numubr is the max texton number +1
histMat=zeros(Nsp,NumOfBins);
    for i = 1:Nsp  
        Indx=seg_k{i};
        histMat(i,:)=hist(texton(Indx),0:NumOfBins-1);
    end
%% compute the chi square distance
total = sum(histMat,2);
total(total==0)=1;
histMat=bsxfun(@rdivide,histMat,total); %normalization
distMat = zeros(Nsp,Nsp);
for i = 1:Nsp
    temp1 = bsxfun(@minus, histMat,histMat(i,:));
    temp2 = bsxfun(@plus, histMat,histMat(i,:));
    temp2(temp2==0)=1; % set those bin = 0 in both two histogram to 1. In this case temp1 also equal to 0, so, the result temp1.^0.5./temp2==0
    temp3 = (temp1.^2)./temp2;
    distMat(i,:) = sum(temp3,2)*0.5;
end

