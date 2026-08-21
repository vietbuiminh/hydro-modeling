using CairoMakie

struct PorosityWindow
    cross::Matrix{Bool}
end

function window_ranges(pw::PorosityWindow, size::Int)
    height, width = Base.size(pw.cross)
    # need to check that size is valid
    # this would be a good place to use a @assert, 
    # but I want to throw an ArgumentError instead
    # comment the checing section to speed up the function
    # 0 < size <= min(height, width) ||
    #     throw(ArgumentError("window size $size must be in 1:$(min(height, width))"))
    half = size ÷ 2
    row_start = (height ÷ 2) - half + 1
    col_start = (width ÷ 2) - half + 1
    return row_start:(row_start + size - 1), col_start:(col_start + size - 1)
end

function preview_window(pw::PorosityWindow, size::Int;
    void_color = :white, solid_color = :black)
    rows, cols = window_ranges(pw, size)
    fig = Figure()
    ax = Axis(
        fig[1, 1],
        title = "$(size) x $(size)",
        aspect = DataAspect(),
    )
    heatmap!(ax, pw.cross, colormap = [void_color, solid_color])
    poly!(
        ax,
        Point2f[
            (rows[1] - 1, cols[1] - 1), (rows[end], cols[1] - 1),
            (rows[end], cols[end]), (rows[1] - 1, cols[end]),
        ],
        color = (:red, 0.0), strokecolor = :red, strokewidth = 3,
    )
    display(fig)
    return fig
end

function calc_avg_porosity(pw::PorosityWindow, size::Int)
    rows, cols = window_ranges(pw, size)
    window = @view pw.cross[rows, cols]
    bulk = length(window)
    void = Int(count(!, window))
    ϕ = void / bulk * 100
    return ϕ
end
