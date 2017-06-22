function F_Dist = compute_FoerstnerDist(Cov_Matrix,epsilon)

if nargin<2
    epsilon=1e-5;
end

K = length(Cov_Matrix);
temp=Cov_Matrix{1};
m=size(temp,1);
F_Dist = zeros(K);

Epsilon=ones(m,1)*epsilon; clear temp m ;
Epsilon=diag(Epsilon);

for i = 1:K    
    Mat_A = Cov_Matrix{i}+Epsilon;
    for j = 1:i
        Mat_B = Cov_Matrix{j}+Epsilon;
        [~,lambda]=eig(Mat_A,Mat_B,'chol');
        lambda=diag(lambda);
        F_Dist(i,j) = sqrt(sum(log(lambda).^2));
    end
  %  F_Dist(i,:) = normalize(F_Dist(i,:)'); % to make the value into [0,1]
end
F_Dist = F_Dist + F_Dist'; % Notice: the diagonal is zero, So, no need to minus the dupilicated diagonal. 
 % to make the value into [0,1]
 %F_Dist = F_Dist/(max(max(F_Dist))-min(min(F_Dist)));
 %F_Dist=(F_Dist+F_dist')*0.5;

end


