function sim_output = simColor(seg_intensity_vals,nb)

NumOfSp = size(seg_intensity_vals,1);

if size(seg_intensity_vals,2)==3
    temp1 = reshape(seg_intensity_vals,[NumOfSp,1,3]);
    temp1 = repmat(temp1,1,NumOfSp);
    temp2 = cat(3,seg_intensity_vals(:,1)',seg_intensity_vals(:,2)',seg_intensity_vals(:,3)');
    temp2 = repmat(temp2,NumOfSp,1);
    dist = sqrt(sum((temp1-temp2).^2,3));
elseif size(seg_intensity_vals,2)==1
    temp1 = repmat(seg_intensity_vals,1,NumOfSp);
    temp2 = repmat(seg_intensity_vals,NumOfSp,1);
    dist = sqrt((temp1-temp2).^2);
end
norm_dist = normalize(reshape(dist,NumOfSp*NumOfSp,1));
norm_dist = reshape(norm_dist, [NumOfSp,NumOfSp]);
sim =exp(-norm_dist); %% use exponetional function convert the distance to similarity,this sim is symetric.
sim_output = sim -diag(diag(sim)); % remove the diagonal, because, in this matrix, the max value is itself.
%sim_output = prune_knn(sim,nb);

%[Y,I] = max(norm_sim,[],2);
%[temp,originalpos] = sort(norm_sim,2,'ascend');
%sim = zeros(NumOfSp,NumOfSp);
% for i = 1:NumOfSp
%     for j = 1:nb    %%find the n closest one; 
%         sim(i,originalpos(i,j)) = temp(i,j)*SizeOfSp(originalpos(i,j));
%     end
% end

end
    


