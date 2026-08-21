using CairoMakie

struct PorosityWindow
    cross::Matrix{Bool}
    label::String
end

function window_ranges(pw::PorosityWindow, width::Int)
    height, w = size(pw.cross)
    0 < width <= min(height, w) ||
        throw(ArgumentError("window width $width must be in 1:$(min(height, w))"))
    half = width ÷ 2
    row_start = (height ÷ 2) - half + 1
    col_start = (w ÷ 2) - half + 1
    return row_start:(row_start + width - 1), col_start:(col_start + width - 1)
end

function preview_window(pw::PorosityWindow, width::Int;
    void_color = colorant"#F0DC82", solid_color = colorant"#70673b")
    rows, cols = window_ranges(pw, width)
    fig = Figure()
    ax = Axis(
        fig[1, 1],
        title = "$(pw.label) — $(width)×$(width) window",
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

function calc_avg_porosity(pw::PorosityWindow, width::Int)
    rows, cols = window_ranges(pw, width)
    window = @view pw.cross[rows, cols]
    return count(!, window) / length(window) * 100
end
