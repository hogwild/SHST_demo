function LogE_Dist = compute_LogEuclideanDist(Cov_Matrix,epsilon)

if nargin<2
    epsilon=1e-5;
end

K = length(Cov_Matrix);
%temp=Cov_Matrix{1};
%m=size(temp,1);
LogE_Dist = zeros(K);
[m,~]=size(Cov_Matrix{1});
%Epsilon=ones(m,1)*epsilon; clear temp m ;
%Epsilon=diag(Epsilon);

for i = 1:K    
    Mat_A = Cov_Matrix{i}+epsilon*eye(m);
    for j=1:i %j = 1:K %because the distance matrix is symetric.It can reduce half of the iterations.
        Mat_B = Cov_Matrix{j}+epsilon*eye(m);
        
        LogE_Dist(i,j) = norm(logm(Mat_A)-logm(Mat_B),'fro'); %Frobenius norm.
    end
    %LogE_Dist(i,:) = normalize(LogE_Dist(i,:)'); % to make the value into [0,1]
    
end
LogE_Dist=LogE_Dist+LogE_Dist'; % Notice: the diagonal is zero, So, no need to minus the dupilicated diagonal. 
 % to make the value into [0,1]
 %F_Dist = F_Dist/(max(max(F_Dist))-min(min(F_Dist)));
 %F_Dist=(F_Dist+F_dist')*0.5;

end


