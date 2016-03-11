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

    # mayavi_plotting.jl
    figure, plot_points, axes, sync,

    # register.jl
    register, iterative_register, register_exact,

    # synchronization.jl
    maximize_trace, parallel_pairwise_register, project_to_stiefel, stiefel_to_permutations,

    # pointsets.jl
    trim, random_perturb!, random_permute!, random_rotate!, pointsets_transform, pointsets_rmsd

include("accelerated_routines.jl")
include("clustering.jl")
include("grid.jl")
include("GPM.jl")
include("density_grid.jl")
include("mayavi_plotting.jl")
include("register.jl")
include("synchronization.jl")
include("pointsets.jl")

end
