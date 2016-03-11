function figure(fig_num::Int, size = 700, fgcolor = (1.0, 1.0, 1.0), bgcolor = (0.0, 0.0, 0.0))
    fig = mlab.figure(
        fig_num,
        size = (size, size),
        fgcolor = fgcolor,
        bgcolor = bgcolor
    )

    mlab.clf(fig)
    return fig
end

function plot_points{T}(fig_num::Int, A::Matrix{T}, color = (0.0, 1.0, 0.0), scale_factor = 0.3, bw = 12.0)
    fig = mlab.figure(fig_num)

    mlab.points3d(
    	figure = fig,
        [0.0], [0.0], [0.0],
        color = (1.0, 0.2, 0.2),
        mode = "sphere",
        scale_factor = 0.2,
	resolution = 12,
        reset_zoom = false
    )

    mlab.points3d(
	figure = fig,
        A[1, :], A[2, :], A[3, :],
        color = color,
        mode = "sphere",
        scale_factor = scale_factor,
	resolution = 12,
        reset_zoom = false
    )
end

function axes(fig_num::Int)
    fig = mlab.figure(fig_num)
    mlab.orientation_axes(figure = fig)
end

function sync(fig_num_ref::Int, fig_num_target::Int)
    mlab.sync_camera(mlab.figure(fig_num_ref), mlab.figure(fig_num_target))
end
