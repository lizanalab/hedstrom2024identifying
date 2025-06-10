include("communities.jl")
using Statistics; mean
using Random; shuffle
using ProgressMeter


# Function to count how many times each community in the original community vector 
# is repeated or matched in the noise communities, based on a similarity threshold.
# Arguments:
#   - communities: a `Communities` object containing original and noise communities.
#   - similarity_cutoff: threshold to decide if two communities are similar.
#   - method: similarity method to be used (:ji, :overlap, :dice).
# Returns:
#   - counts: a vector where each entry represents the number of times a community is similar across the noise communities.
function count_all_community_repeats(communities::Communities, similarity_cutoff::Float64, method=:ji)
    if method ∉ [:ji, :overlap, :dice, :gamma, :AR]
        throw(ValueError(method, "not a defined method"))
    end

    # Returns a vector that gives the community stability per community. I.e. index is community number
    n_samps = length(communities.noise_comms)
    counts = zeros(length(unique(communities.orig_comm[communities.orig_comm .> 0])))  # Only consider non-zero communities.
    
    # Loop through each noise community sample and compute the max similarity for each community.
    for samp_i in 1:n_samps
        max_similarities = get_max_similarity_between_communities(communities.orig_comm, communities.noise_comms[samp_i], method)
        counts .+= (max_similarities .>= similarity_cutoff)  # Count similarities that exceed the cutoff.
    end
    return counts
end


# Function to calculate the maximum similarity between each original community and noise community.
# Arguments:
#   - orig_comm: vector of the original community assignments.
#   - noise_comm: vector of the noise community assignments.
#   - method: similarity metric (:ji, :overlap, :dice).
# Returns:
#   - max_similarities: a vector of maximum similarities for each original community.
function get_max_similarity_between_communities(orig_comm::Vector{Int}, noise_comm::Vector{Int}, method)
    max_similarities = zeros(length(unique(orig_comm[orig_comm .> 0])))  # Initialize vector to store max similarities.
    used_comm_j = zeros(length(max_similarities))  # Track which noise communities are already assigned.

    # Loop through each unique original community and compare it to noise communities.
    similarity = 0
    for comm_i in unique(orig_comm[orig_comm .> 0]) 
            # Compute similarity based on the chosen method.
        for comm_j in unique(noise_comm[orig_comm .== comm_i])
            w11 = sum(orig_comm .== comm_i .&& noise_comm .== comm_j)
            w10 = sum(orig_comm .== comm_i .&& noise_comm .!== comm_j)
            w01 = sum(orig_comm .!== comm_i .&& noise_comm .== comm_j)
            w00 = sum(orig_comm .!== comm_i .&& noise_comm .!== comm_j)
            M = w11 + w01 + w10 + w00
            if method == :ji
                # similarity = sum(orig_comm .== comm_i .&& noise_comm .== comm_j)/sum(orig_comm .== comm_i .|| noise_comm .== comm_j)
                similarity = w11/(w10+w01+w11)
            elseif method == :overlap
                similarity = sum(orig_comm .== comm_i .&& noise_comm .== comm_j)/min(sum(orig_comm .== comm_i), sum(noise_comm .== comm_j))
            elseif method == :dice
                similarity = 2*sum(orig_comm .== comm_i .&& noise_comm .== comm_j)/(sum(orig_comm .== comm_i) + sum(noise_comm .== comm_j))
            elseif method == :gamma
                similarity = (M*w11 - (w11 + w10)*(w11 + w01))/(sqrt((w11 + w10)*(w11 + w01)*(M - (w11 + w10))*(M - (w11 + w01))))
            elseif method == :AR
                similarity = (w11 - 1/M*(w11 + w10)*(w11 + w01))/(0.5*((w11 + w10) + (w11 + w01)) - 1/M*(w11 + w10)*(w11 + w01))
            else
                throw(ValueError(method, "not a defined method"))
            end
            # Update the max similarity if a better match is found and the community isn't used yet.
            if similarity > max_similarities[comm_i] && (comm_j ∉ used_comm_j)
                max_similarities[comm_i] = similarity
                used_comm_j[comm_i] = comm_j  # Mark the noise community as used.
            end
        end
    end
    return max_similarities  # Return the vector of maximum similarities for each original community.
end


# Function to compute the percentage overlap of community borders between different noise communities.
# The overlap is scaled by the number of borders compared to the original community border set.
# Arguments:
#   - communities: a `Communities` object containing original and noise communities.
# Returns:
#   - percentage_overlap: the average percentage overlap between community borders across all noise samples.
function community_border_overlap(communities::Communities)
    # Calculates the overlap (JI) between all community borders scaled by the number of borders compared to original border set
    borders_orig = find_community_borders(communities.orig_comm)
    n_samps = length(communities.noise_comms)  # Number of noise community samples.
    percentage_overlap = 0.0
    # Loop through pairs of noise community samples and compute border overlaps.
    for comm_i in 1:n_samps
        borders_i = find_community_borders(communities.noise_comms[comm_i])  # Borders in the current noise community.
        for comm_j in comm_i+1:n_samps
            borders_j = find_community_borders(communities.noise_comms[comm_j])  # Borders in another noise community.

            # Calculate union and intersection of borders and the scaled overlap.
            border_union = length(union(borders_i, borders_j))
            border_min = min(length(borders_i), length(borders_j))
            percentage_overlap += length(intersect(borders_i, borders_j))/border_union * (1 - abs(border_min - length(borders_orig))/border_min)
        end
    end
    # Normalize the overlap by the number of pairwise comparisons.
    percentage_overlap /= (n_samps * (n_samps - 1) / 2)
    return percentage_overlap
end


# Function to find the indices where the community changes, i.e., the borders between communities.
# Arguments:
#   - communities: vector of community assignments.
# Returns:
#   - A vector of indices where community transitions occur.
function find_community_borders(communities::Vector{Int})
    # Find where the community assignments change between consecutive elements and return those indices.
    return findall(diff(communities) .!= 0) .+ 1
end
