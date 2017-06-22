function sim_output = simTextons(seg,nb,textonPath,textonDistPath,k)
NumOfSp = length(seg);
idxOfTextonDistMat = int2str(k);
idxOfTextonDistMat2 = int2str(NumOfSp);
thePath = [textonDistPath,'\textonDist',idxOfTextonDistMat,'_NumOfSp',idxOfTextonDistMat2,'.mat'];
if ~exist(thePath,'file')
    if~exist(textonDistPath,'dir')
        mkdir(textonDistPath);
    end
    if ~exist(textonPath,'file')
        assert(false,'Error: Please select a distance for covariance matrix.');
    else
        load(textonPath); % load the textion matrix
        distMat = compute_ChiSqDist(texton,seg);
        %save distance
        save(thePath,'distMat');
    end
else
%load distance
    load(thePath);
end

% note the distance to the superpixel itself is always zero, so the similarest one is
% itself, and the diagonal of distMat is 0.
norm_dist = normalize(distMat(:)); % the distance is normalized on all entries in the matrix
norm_dist = reshape(norm_dist,[NumOfSp,NumOfSp]);
sim = exp(-norm_dist);
sim_output = sim-diag(diag(sim)); %remove the diagonal set the similarity to itself to be zero.
% just keep the most similariy one
%sim_output = prune_knn(sim,nb); 