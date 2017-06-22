function sim_output = simCovMatrix(img,seg,DistType,nb,covDispath,k)
NumOfSp = length(seg);
idxOfCovMat = int2str(k);
idxOfCovMat2 = int2str(NumOfSp);

if strcmp(DistType,'LogEuclidean')
    thePath = [covDispath,'\LogEuclidean',idxOfCovMat,'_NumOfSp',idxOfCovMat2,'.mat'];
    if ~exist(thePath,'file')
        if~exist(covDispath,'dir')
            mkdir(covDispath);
        end
        covMat = compute_covMatrix(img,seg); %%compute the covariance Matrix   
        distMat = compute_LogEuclideanDist(covMat);
    %save distance
        save(thePath,'distMat');
    else
    %load distance
        load(thePath);
    end
elseif strcmp(DistType,'Foerstner')
    thePath = [covDispath,'\Foerstner',idxOfCovMat,'_NumOfSp',idxOfCovMat2,'.mat'];
    if ~exist(thePath,'file')
        if~exist(covDispath,'dir')
            mkdir(covDispath);
        end
        covMat = compute_covMatrix(img,seg); %%compute the covariance Matrix   
        distMat = compute_FoerstnerDist(covMat);
    %save distance
        save(thePath,'distMat');
    else
    %load distance
        load(thePath);
    end
else
    assert(false,'Error: Please select a distance for covariance matrix.');
end

% note the distance to the superpixel itself is always zero, so the similarest one is
% itself, and the diagonal of distMat is 0.
norm_dist = normalize(distMat(:)); % the distance is normalized on all entries in the matrix
norm_dist = reshape(norm_dist,[NumOfSp,NumOfSp]);
sim = exp(-norm_dist);
sim_output = sim-diag(diag(sim)); %remove the diagonal set the similarity to itself to be zero.
% just keep the most similariy one
%sim_output = prune_knn(sim,nb); 

