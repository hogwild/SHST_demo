function sim_output=simColorSAS(vals,edges,valScale,EPSILON)

%Constants
if nargin < 4
    EPSILON = 1e-5;
end

%Compute intensity differences
if valScale > 0
    valDistances=sqrt(sum((vals(edges(:,1),:)- ...
        vals(edges(:,2),:)).^2,2));
    valDistances=normalize(valDistances); %Normalize to [0,1]
else
    valDistances=zeros(size(edges,1),1);
    valScale=0;
end
%Compute the neighbourhood matrix
sim_output=zeros(size(vals,1));

%Compute Gaussian weights
sim=exp(-valScale*valDistances)+...
     EPSILON;
 % convert edges to neighbourhood in matrix
 ind = sub2ind([size(vals,1),size(vals,1)],[edges(:,1);edges(:,2)],[edges(:,2);edges(:,1)]);
 sim_output(ind) = [sim;sim];
 %sim_output = sim_output.*RegionSize;

% for i=1:size(edges(:,1))
%     sim_output(edges(i,1),edges(i,2))=sim(i);
%     sim_output(edges(i,2),edges(i,1))=sim(i);
% end

