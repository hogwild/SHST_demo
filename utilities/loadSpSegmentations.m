%% some Auxilary functions
% load superpixel segmentations (for using parfor)
function [seg,labels_img,seg_vals,seg_lab_vals,seg_edges,seg_img]=loadSpSegmentations(SpPath)
load(SpPath);
end