function SimMat = SpAssoc_SimMeasure(HistArray,RegionSize)
%% this function compute the similarity between superpixel associations
% Input:
% HistArray - an array of histogram representation of the superpixel
% associations,each row represent the histogram of one sp assoc.
% RegionSize -  the size of the superpixe associations, it is a vector
% Output:
% sim - the similarity matrix of Sp Assoc
%% main
NumOfSpAssoc = length(RegionSize);
SimMat = zeros(NumOfSpAssoc);
for i = 1:NumOfSpAssoc
    HistInnerProduct = HistArray*HistArray(i,:)';
    Size = RegionSize*RegionSize(i);
    SimMat(i,:) = HistInnerProduct./Size;
end
% HistInnerProduct = HistArray*currentHist;
% Size = RegionSize*CurrentRegionSize;
% SimMat = HistInnerProduct./Size;
