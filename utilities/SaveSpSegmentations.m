%% some Auxilary functions
% save superpixel segmentations (for using parfor)
function SaveSpSegmentations(SaveSpPath,seg,labels_img,seg_vals,seg_lab_vals,seg_edges,seg_img)
% Seg = seg;
% Labels_img = labels_img;
% Seg_vals = seg_vals;
% Seg_lab_vals = seg_lab_vals;
% Seg_edges = seg_edges;
% Seg_img = seg_img;
save(SaveSpPath, 'seg', 'labels_img', 'seg_vals', 'seg_lab_vals' ,'seg_edges','seg_img');
end