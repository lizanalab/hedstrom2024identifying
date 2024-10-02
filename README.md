# Stable community detection and analysis

This is the code repository for the paper "Identifying stable communities in Hi-C data using a multifractal null model".

```
@article{hedstrom2024identifying,
  title={Identifying stable communities in Hi-C data using a multifractal null model},
  author={Hedstr{\"o}m, Lucas and Mart{\'\i}nez, Ant{\'o}n Carcedo and Lizana, Ludvig},
  journal={arXiv preprint arXiv:2405.05425},
  year={2024}
}
```

This project provides functions for analyzing communities in a given dataset, with a specific focus on Hi-C data and the impact of noise on community structures. It includes tools for calculating community stability, community border overlaps, and other related operations. The project uses Julia and MATLAB as the primary languages.

This project also uses the MATLAB implementation of the [generalized Louvain method](https://github.com/GenLouvain/GenLouvain);

```
Lucas G. S. Jeub, Marya Bazzi, Inderjit S. Jutla, and Peter J. Mucha, "A generalized Louvain method for community detection implemented in MATLAB," https://github.com/GenLouvain/GenLouvain (2011-2019).
```

## Example data

This repository contains the `.RAWobserved` Hi-C contacts of human chromosome 10 at 100kb resolution. This data is taken from the [Dec. 2013 assembly of the human genome (hg38, GRCh38 Genome Reference Consortium Human Reference 38 (GCA_000001405.15))](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/).

## Usage

### Generating data

The primary function in this project is `generate_communities`, which detects communities from Hi-C data. To use it, run the following in MATLAB:

```matlab
generate_communities(chr_num, resolution, gamma, delta, a, num_samp)
```

#### Parameters:
- `chr_num`: Chromosome number (e.g., `1`, `2`).
- `resolution`: Hi-C data resolution (in base pairs, e.g., `100000` for 100kb).
- `gamma`: Resolution parameter for the community detection algorithm.
- `delta`: Scaling factor for the diagonal noise deviations in the Hi-C matrix.
- `a`: "Structure" parameter for hierarchical domain model.
- `num_samp`: Number of noisy community samples to generate.

#### Example:
```matlab
generate_communities(1, 100000, 1.0, 0.5, 0.3, 100)
```

This will generate communities for chromosome 1, at 100kb resolution, using a gamma of `1.0`, delta of `0.5`, and parameter `a` of `0.3`, producing 100 noisy samples.

#### Input Files

- **Hi-C Raw Data:** The function expects the raw Hi-C contact data in a specific file structure, located at:
  ```
  data/chr{chr_num}_{resolution}kb.RAWobserved
  ```
  - Example: `data/chr10_100kb.RAWobserved`

#### Output Files

The function will generate the following output files:

- `chr_{chr_num}.data.info`: Metadata file with information about the chromosome and matrix dimensions.
- `chr{chr_num}_{resolution}kb_gamma{gamma}_delta{delta}_a{a}.csv`: A CSV file that contains the community assignments for the original Hi-C data and for each noisy sample.

#### CSV Format:

The CSV file contains the community assignments for each sample. The first column is the sample number (`0` for the original data, followed by the noisy samples `1` through `num_samp`), and the subsequent columns represent the community assignment for each node (genomic region).

### Analysing data

The code to analyse data utilises Julia. Below follows a usage example and a description of the different function. Refer to the paper for more information.

#### Example
```julia
# Load necessary functions
include("community_analysis.jl")
include("parse.jl")

# Read Hi-C and community files
hic_matrix = read_hic_file("data/chr1_100kb.RAWobserved", 1000, 100000)
communities = read_community_file("data/community.csv")

# Calculate the overlap of community borders across noise samples
overlap = community_border_overlap(communities)

# Count how many times each community in the original structure repeats in the noise communities
repeat_counts = count_all_community_repeats(communities, 0.5, :ji)
```

### Communities data structure
This custom data structure represents the original community structure and a set of noise community structures. It has the following components:
- `orig_comm`: A vector representing the original community assignments.
- `noise_comms`: A list of noise community vectors (perturbed versions of `orig_comm`).
- Other parameters such as chromosome number, resolution, and clustering parameters.

### Data analysis function description

#### 1. Community Repeat Counting

##### `count_all_community_repeats(communities::Communities, similarity_cutoff::Float64, method=:ji)`

This function calculates how many times each original community is repeated across noise samples. You can specify different methods for measuring similarity, including:
- Jaccard Index (`:ji`)
- Overlap coefficient (`:overlap`)
- Dice coefficient (`:dice`)

The function returns a vector where each entry represents the count of how often the corresponding community in the original structure is repeated across noise samples.

#### 2. Community Similarity Calculation

##### `get_max_similarity_between_communities(orig_comm::Vector{Int}, noise_comm::Vector{Int}, method)`

This function calculates the maximum similarity between each original community and noise communities. It supports multiple similarity measures, such as Jaccard, overlap, and Dice coefficients.

#### 3. Community Border Overlap

##### `community_border_overlap(communities::Communities)`

This function computes the percentage overlap of community borders across all noise community samples. The overlap is scaled by the number of borders compared to the original community borders.

The output is a single percentage representing how similar the borders are across different noise samples.

##### `find_community_borders(communities::Vector{Int})`

This helper function identifies the borders (indices) where communities change in a given vector of community assignments. It is used internally by the `community_border_overlap` function.

#### 4. Centromere Handling

##### `remove_centromere_from_communities!(communities::Communities, hic::Matrix{Float64})`

This function removes the centromeres (non-informative regions) from the community data. Centromeres are identified from the Hi-C contact matrix (`hic`), and their corresponding regions are removed from both the original and noise communities.

##### `get_indices_of_centromeres(hic::Matrix{Float64})`

This function identifies the indices of centromeres in the Hi-C matrix by finding rows and columns with zero interaction sums. These indices are used by the centromere removal functions.

#### 5. Hi-C File Processing

##### `read_hic_file(file_path::AbstractString, num_nodes::Int, res::Int)::Matrix{Float64}`

This function reads a Hi-C interaction file (in a tab-delimited edge list format) and creates a Hi-C matrix. It takes the number of nodes and the resolution as arguments and returns a matrix of contact frequencies between genomic regions.

#### 6. Community File Processing

##### `read_community_file(file_path::AbstractString; remove_centromere_communities=true)::Communities`

This function reads a community file (in a tab-delimited format) and constructs a `Communities` object. The file is expected to follow a specific naming convention to extract parameters such as chromosome number and resolution from the filename. See the output from the MATLAB function.

It optionally removes centromere regions from the communities, using an external Hi-C file for centromere identification. However, this is hard-coded and should be modified for the specific use case.

#### 7. HMM File Processing

##### `read_hmm_file(file_path::AbstractString, chr::Int, num_nodes::Int, res::Int)`

This function reads a Hidden Markov Model (HMM) file and processes state information for a given chromosome. It organizes data based on state-specific information and computes the fold enrichment for each state.

The function returns the processed data, sorted keys representing the HMM states, and fold enrichment values.