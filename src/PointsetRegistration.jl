module PointsetRegistration

using MoreMatrices
using Munkres

export
    # exact_pairwise.jl
    register,

    # permutation_synchronization.jl
    synchronize_permutations

include("exact_pairwise.jl")
include("exact_pairwise_multi.jl")
include("permutation_synchronization.jl")

end
