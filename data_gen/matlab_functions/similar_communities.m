function [similarCommunities,num_similar] = similar_communities(jaccard_matrix,similarity_threshold)
        
    similarCommunities = cell(1, size(jaccard_matrix, 1));
    num_similar = 0;
    for i = 1:size(jaccard_matrix, 1)
        similarCommunities{i} = find(jaccard_matrix(i, :) >= similarity_threshold);
        num_similar = num_similar + numel(similarCommunities{i});
    end

end