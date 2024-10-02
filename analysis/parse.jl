using DelimitedFiles
include("communities.jl")

output_folder = "code/data_gen/output"
data_folder = "data"


# Function to remove centromere regions from communities based on a Hi-C matrix.
# Arguments:
#   - communities: a `Communities` object containing original and noise communities.
#   - hic: a matrix containing Hi-C contact data.
# The function modifies the `communities` object by removing centromere regions.
function remove_centromere_from_communities!(communities::Communities, hic::Matrix{Float64})
    centromere_inds = get_indices_of_centromeres(hic)  # Get centromere indices from the Hi-C matrix.
    remove_centromere_from_community!(communities.orig_comm, centromere_inds)  # Remove from original community.
    for comm_i in eachindex(communities.noise_comms)
        remove_centromere_from_community!(communities.noise_comms[comm_i], centromere_inds)  # Remove from noise communities.
    end
end


# Helper function to remove centromere regions from a specific community.
# Arguments:
#   - comm: a vector representing a community.
#   - centromere_inds: vector of centromere indices to be removed from the community.
# The function modifies the `comm` vector by setting centromere indices to zero.
function remove_centromere_from_community!(comm::Vector{Int}, centromere_inds::Vector{Int})
    map(centromere_inds) do centromere_ind
        zero_val = comm[centromere_ind]  # Get the community value at the centromere index.
        
        # Set centromere community to zero and decrement other community indices if necessary.
        if zero_val !== 0
            comm[centromere_ind] = 0
            comm[comm .> zero_val] .-= 1  # Decrement community indices greater than centromere's value.
        end
    end 
end


# Function to obtain the indices of centromeres from a Hi-C matrix.
# Arguments:
#   - hic: a Hi-C contact matrix.
# Returns:
#   - result_indices: a sorted vector of indices where the row or column sum is zero (indicative of centromere positions).
function get_indices_of_centromeres(hic::Matrix{Float64})
    row_sums = sum(hic, dims=2)  # Sum the values of each row.
    col_sums = sum(hic, dims=1)  # Sum the values of each column.

    # Find the indices where the combined row and column sums equal zero (centromeres).
    result_indices = sort(findall(x -> x == 0, (row_sums' .+ col_sums)[:]))

    return result_indices
end


# Function to read a Hi-C file and return a Hi-C matrix.
# Arguments:
#   - file_path: the path to the file containing the Hi-C data.
#   - num_nodes: the expected number of nodes in the Hi-C matrix.
#   - res: the resolution of the Hi-C data.
# Returns:
#   - hic: a matrix representing the Hi-C contact data.
function read_hic_file(file_path::AbstractString, num_nodes, res)::Matrix{Float64}
    # Read the edge list (source, target, weight) from the file.
    edge_list = readdlm(file_path, '\t', Float64, '\n')

    # Initialize an empty Hi-C matrix with zeros.
    hic = zeros(Float64, num_nodes, num_nodes)

    max_ind = 1
    # Populate the Hi-C matrix with the edge weights.
    for i in 1:size(edge_list)[1]
        edge = edge_list[i, :]
        source, target, weight = edge
        x = Int(source/res) + 1
        y = Int(target/res) + 1
        hic[x, y] = Float64(weight)  # Set both symmetric positions in the matrix.
        hic[y, x] = Float64(weight)
        
        # Track the maximum index used to resize the matrix if needed.
        if x > max_ind
            max_ind = x
        end
        if y > max_ind
            max_ind = y
        end
    end

    # If the actual matrix size is smaller than expected, truncate it to the correct size.
    if max_ind < num_nodes
        hic = hic[1:max_ind, 1:max_ind]
    end

    return hic
end


# Function to read a community file and construct a `Communities` object.
# Arguments:
#   - file_path: path to the file containing community data.
#   - remove_centromere_communities: boolean flag to indicate whether centromere communities should be removed (default is true).
# Returns:
#   - communities: a `Communities` object containing the community data.
function read_community_file(file_path::AbstractString; remove_centromere_communities=true)::Communities
    # Extract metadata (chromosome number, resolution, gamma, delta, and a) from the file name.
    file_name_regex = match(r"chr(\d+)_(\d+)kb_gamma([\d.]+)_delta([\d.]+)_a([\d.]+)\.csv", basename(file_path))

    # Read the community data from the file.
    communities = readdlm(file_path, '\t', Float64, '\n', skipstart=1)
    
    N = length(communities[2, 2:end])  # Number of nodes.
    chr = parse(Int, file_name_regex[1])  # Chromosome number.
    res = parse(Int, file_name_regex[2])*1000  # Resolution.

    # Construct the `Communities` object with metadata and the community data.
    communities = Communities(
        N,
        chr,
        res,
        parse(Float64, file_name_regex[3]),  # Gamma.
        parse(Float64, file_name_regex[4]),  # Delta.
        parse(Float64, file_name_regex[5]),  # a parameter.
        communities[1, 2:end],  # Original community data.
        Dict(i => communities[i+1, 2:end] for i in 1:size(communities)[1]-1)  # Noise communities.
    )

    # If centromere communities should be removed, process the communities.
    if remove_centromere_communities
        hic = read_hic_file("$(data_folder)/chr$(chr)_100kb.RAWobserved", N, res)
        remove_centromere_from_communities!(communities, hic)
    end

    return communities
end


# Function to read an HMM (Hidden Markov Model) file of chromatin states and compute data for peak analysis.
# Arguments:
#   - file_path: path to the HMM file.
#   - chr: chromosome number to filter data.
#   - num_nodes: the number of nodes (bins) in the Hi-C matrix.
#   - res: resolution of the data (bin size).
# Returns:
#   - data: a dictionary where keys are states and values are the associated data vectors.
#   - sorted_keys: a list of sorted state keys.
#   - fex: a dictionary of fold enrichment data.
function read_hmm_file(file_path::AbstractString, chr, num_nodes, res)
    hmm = readdlm(file_path, '\t', String, '\n', skipstart=1)  # Read the HMM data.
    data = Dict()  # Dictionary to store data per state.
    kx = Dict()  # Dictionary to store the number of peaks per bin.

    N = size(hmm)[1] / length(unique(hmm[:, 2]))  # Calculate total number of peaks per chromosome.
    
    # Process each row in the HMM file.
    for row_i in 1:size(hmm)[1]
        hmm_row = hmm[row_i, :]
        chrom = hmm_row[2]  # Chromosome identifier.

        # Skip rows not corresponding to the specified chromosome.
        (chrom !== "chr$(chr)") && continue

        state = hmm_row[5]  # State identifier.
        start = parse(Int, hmm_row[3])  # Start position.
        stop = parse(Int, hmm_row[4])  # Stop position.

        # Allocate data and peak vectors if they do not exist for the state.
        (!haskey(data, state)) && (data[state] = zeros(num_nodes))
        (!haskey(kx, state)) && (kx[state] = zeros(num_nodes))

        bin_start = Int((start - start % res) / res) + 1  # Calculate bin index.

        # Increment the peak count for the bin if it is within bounds.
        if bin_start < num_nodes
            kx[state][bin_start] += 1
        end

        # Distribute the data across the bins until the stop position is reached.
        while start < stop
            nextend = min(bin_start * res, stop)
            if bin_start > num_nodes
                println("for chr$(chr) we got $bin_start > $num_nodes")
                break
            end
            data[state][bin_start] += nextend - start
            start = bin_start * res
            bin_start += 1
        end
    end

    # Calculate fold enrichment per state.
    n = reduce(.+, values(kx))
    fex = Dict()
    for (key, vector) in kx
        KX = sum(vector)
        fex[key] = vector ./ (KX .* (n / N))  # Compute fold enrichment.
        fex[key][isnan.(fex[key])] .= 0  # Set NaNs to 0.
    end

    # Normalize data by resolution.
    for (key, vector) in data
        data[key] = vector ./ res
    end

    # Custom comparison function to sort the keys by state number.
    function hmm_comp(key1, key2)
        return parse(Int, split(key1, "_")[1]) < parse(Int, split(key2, "_")[1])
    end

    sorted_keys = sort(collect(keys(data)), lt=hmm_comp)  # Sort the state keys.

    return data, sorted_keys, fex  # Return the processed data, sorted keys, and fold enrichment.
end
