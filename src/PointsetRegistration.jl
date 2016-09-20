module PointsetRegistration

export
    # accelerated_routines.jl
    neg_XtRY, wYPXt, cost, partial_cost,

    # clustering.jl
    kernel_kmeans, max_k_cut, factor_k_cut,

    # grid.jl
    grid, grid_coords,

    # GPM.jl
    GPM_orthogonal, GPM_permutation,

    # density_grid.jl
    parallel_density_grid, density_grid, density_grid_stats,

    # register.jl
    register, exact_register, gaussian_register, register_to_reference,

    # synchronization.jl
    parallel_pairwise_register, stiefel_to_permutations,

    # pointsets.jl
    trim, random_perturb!, random_permute!, random_rotate!, random_rotate_reflect!, pointsets_transform, pointsets_rmsd,

    # group_representations
    irrep_O3!, irrep_O3


include("accelerated_routines.jl")
include("clustering.jl")
include("grid.jl")
include("GPM.jl")
include("density_grid.jl")
include("register.jl")
include("synchronization.jl")
include("pointsets.jl")
include("irreps_O3.jl")

end
