function [NumOfTrees,rootId] = Check_Graph_Connection(Graph)
%% this function find the connection of a graph by Warshell algorithm
% Input:
% Graph - the graph need to be check
% Output:
% Connection - the number of trees  of the graph,
C = Graph~=0;
n = size(C,1);
NumOfTrees = 0;
Idx =1;
TreeId = zeros(n,1);
for i = 1:n
    for j = (i+1):n
        if C(i,j)==1
            if TreeId(i) == TreeId(j)
                if TreeId(i)==0
                    TreeId(i) = Idx;
                    TreeId(j) = Idx;
                    Idx = Idx+1;
                    NumOfTrees=NumOfTrees+1;
                end
            else
                if TreeId(i)==0
                    TreeId(i) = TreeId(j);
                elseif TreeId(j)==0
                    TreeId(j) = TreeId(i);
                else
                    TreeId(TreeId==TreeId(i))=TreeId(j); % if the chain with different label are connected, merge them
                    NumOfTrees = NumOfTrees-1;
                end
            end
        end
    end
end
rootId = zeros(NumOfTrees,1);
for i = 1:NumOfTrees
    temp = find(TreeId == i);
    rootId(i) = temp(1);
end

                    
                    
        

