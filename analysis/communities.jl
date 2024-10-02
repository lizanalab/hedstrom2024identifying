struct Communities
    N::Int
    chr::Int
    res::Int
    gamma::Float64
    delta::Float64
    a::Float64
    orig_comm::Vector{Int}
    noise_comms::Dict{Int, Vector{Int}}
end