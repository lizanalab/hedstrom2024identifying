function jaccard_indices = jaccard_communities(communityAssignment1,communityAssignment2)

    numNodes = length(communityAssignment1);
    
    partition1 = cell(max(communityAssignment1), 1);
    partition2 = cell(max(communityAssignment2), 1);

    for node = 1:numNodes
        community1 = communityAssignment1(node);
        partition1{community1} = [partition1{community1}, node];
    
        community2 = communityAssignment2(node);
        partition2{community2} = [partition2{community2}, node];
    end

    numCommunities1 = numel(partition1);
    numCommunities2 = numel(partition2);
   
    jaccard_indices = zeros(numCommunities1, numCommunities2);
    
    for i = 1:numCommunities1
        for j = 1:numCommunities2
            intersection_comm = numel(intersect(partition1{i}, partition2{j}));
            union_comm = numel(union(partition1{i}, partition2{j}));
            jaccard_index = intersection_comm / union_comm;
            jaccard_indices(i,j) = jaccard_index;
        end
    end
end