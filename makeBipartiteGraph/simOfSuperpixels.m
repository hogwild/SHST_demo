function sim = simOfSuperpixels(img,seg, seg_lab_vals, seg_edges, para,featPath,k)
%% This function compute the similarity between the given superpixels 
% Input:
% seg -- the set of oversegmentation labels, i.e. superpixel labels of each pixel in each segmentation;
% seg_lab_val -- the Lab Color value of each superpixel;
% seg_edges -- the adjacency of the superpixels;
% para -- the parameters;
% Output:
% sim -- a cell of (n x n) symetric matrix, the non-zero entris represent the
% weights on the hyper edges.
%%
switch para.mode
    case 'SAS'
        %measure the similarity between superpixels in the same way of
        %SAS, i.e. the color different between adjoint superpixels only
        sim = simColorSAS(seg_lab_vals,seg_edges,para.beta); % in order to compare the efficiency of Vcut       
    case 'ColorCovMat' %the method used in ICIP2014 paper      
        %measure the similarity between superpixels in the same way of
        %Color and Covariance Matrix
        sim_color = simColor(seg_lab_vals,para.nb);  
        sim_texture = simCovMatrix(img,seg,para.distType,para.nb,featPath.CovMat,k); % parameter k is for load the precomputed covmat distance
        sim = 0.5*(sim_color + sim_texture);        
    case 'ColorTexton'
        sim_color = simColor(seg_lab_vals);
        sim_texture = simTextons(seg,para.nb,featPath.Texton,featPath.TextonDist,k);
        sim = 0.5*(sim_color + sim_texture);       
    case 'Neighbours'
        rowsub = [seg_edges(:,1);seg_edges(:,2)];
        colsub = [seg_edges(:,2);seg_edges(:,1)];
        Nsp = size(seg_lab_vals,1);
        ind = sub2ind([Nsp,Nsp],rowsub,colsub);
        neighourMat = zeros(Nsp,Nsp);
        neighourMat(ind) = 1; % the neighourhood relation matrix. For those adjacency superpixel (i,j) is 1, others are 0;
        sim_color = simColor(seg_lab_vals,para.nb);
        sim_texture = simTextons(seg,para.nb,featPath.Texton,featPath.TextonDist,k);
        sim = 0.5*(sim_color + sim_texture).*neighourMat; % pick out the adjacent superpixels
end
