#=
    Author: Viet M. Bui
    Date: 2026-08-21
    GLY 6826 - Hydrogeologic Modeling
    Homework 1: Porosity and Representative Elementary Volume (REV)

    Note: in order to fun the file, make sure to check for `data` 
    and `funcs` folders in the same dir
=#
include("funcs/read_bn_dat.jl")
include("funcs/calc_avg_porosity.jl") # this is an Object
using CairoMakie
using Colors

light_buff = colorant"#F0DC82"
buff = colorant"#70673b"
top_color = colorant"#2c489c"
bot_color = colorant"#a8761e"

set_theme!(theme_latexfonts())
update_theme!(
    fontsize = 20,
    figure_padding = (10, 50, 10, 10),
)

#=
The files berea_xsection_top.dat and berea_xsection_bot.dat
each contain a 400 x 400 pixel matrix representing a cross section 
of a micro-CT image of Berea sandstone. A value of 0 represents 
pore space (void), and a value of 1 represents solid rock.
=#
top_cross =
    read_bn_dat("data/berea_xsection_top.dat", width = 400, height = 400, T =
    Bool) # read in top
bot_cross =
    read_bn_dat("data/berea_xsection_bot.dat", width = 400, height = 400, T =
    Bool) # read in bot
#= expected out:
┌ Info: data/berea_xsection_top.dat:
│   void = 32213
│   void_percentage = "20.13%"
│   solid = 127787
└   solid_percentage = "79.87%"
┌ Info: data/berea_xsection_bot.dat:
│   void = 25343
│   void_percentage = "15.84%"
│   solid = 134657
└   solid_percentage = "84.16%"
=#

#= 
    a) Plot the two micro-CT cross sections.
=#
sandx = Figure(size = (500, 1000), dpi = 300, title = "Berea Sandstone XSection")
ax_top = Axis(sandx[2, 1], title = "Top")
ax_bot = Axis(sandx[3, 1], title = "Bottom")

Label(
    sandx[0, :], text = "Berea Sandstone XSection",
    font = :bold,
    fontsize = 23, tellwidth = false,
)
Legend(
    sandx[1, :],
    [PolyElement(color = light_buff), PolyElement(color = buff)],
    ["Void", "Solid"],
    orientation = :horizontal, tellwidth = false,
)

heatmap!(ax_top, top_cross .>= true, colormap = [light_buff, buff])
heatmap!(ax_bot, bot_cross .>= true, colormap = [light_buff, buff])

save("figures/berea_sandstone_xsection.png", sandx)
display(sandx) # view the plot

#=
    b) Estimate the porosity of the rock from each cross section.
    some of the info has been precal inside the func read_bn_dat,
    but I can also calculate it here for the hw.
    I can use pre cals for checking. 
=#
total_volume = 400*400 #square pixels
void_volume_top = length(findall(x -> x == false, top_cross))
ϕ_top = void_volume_top / total_volume * 100
void_volume_bot = length(findall(x -> x == false, bot_cross))
ϕ_bot = void_volume_bot / total_volume * 100

@info(
    "Estimate ϕ of the rock for each cross section:",
    ϕ_Top = ϕ_top,
    ϕ_Bottom = ϕ_bot
)

#=
    c) For each cross section,
    calculate the average porosity within a square window centered 
    on the image as a function of window size. Plot the results for 
    both cross sections on the same graph.

    Note: I made into an Object/Module so it is easier to debug 
    and reuse.
=#

top_cross_window = PorosityWindow(top_cross)
bot_cross_window = PorosityWindow(bot_cross)

# preview_window(top_cross_window, 100) # uncomment to test a window

widths = 2:2:400 # skip by 2 because window size should be even
ϕ_top_by_width = [calc_avg_porosity(top_cross_window, w) for w in widths]
ϕ_bot_by_width = [calc_avg_porosity(bot_cross_window, w) for w in widths]

rev = Figure(size = (700, 500))
ax_rev = Axis(
    rev[1, 1],
    title = "ϕ as a function of window size",
    xlabel = "Window size (pixels)",
    ylabel = "Average ϕ (%)",
)
lines!(ax_rev, widths, ϕ_top_by_width, label = "Top", linewidth = 2, color = top_color)
scatter!(ax_rev, widths, ϕ_top_by_width, color = top_color)
lines!(
    ax_rev,
    widths,
    ϕ_bot_by_width,
    label = "Bottom",
    linewidth = 2,
    color = bot_color,
)
scatter!(ax_rev, widths, ϕ_bot_by_width, color = bot_color)
rev_legend = axislegend(ax_rev, position = :cb)
save("figures/berea_sandstone_avg_phi.png", rev)
display(rev) # view the plot

#=
    d) Based on the results from part (c), estimate the characteristic length
    scale of the representative elementary volume (REV) in pixel units.

    Note: the REV is the large enough window size which the average ϕ 
    properties become approx. stable. Also it is better to let the computer 
    do this by scanning the part (c) curve from the largest to smallest 
    windows and take the first window size after the last tolerance violation.

    After running, REV estimated for Top xcross is 344 px (ϕ ~= 20%) because it
    reached a stable flat line in the range of 344 to 400 px,
    and for Bot xcross is 300 px (ϕ ~= 16%). However, the Bot curve suggest the 
    REV could be improved by even a larger window size since the curve is still
    going down and have not yet reached a stable flat line. 
=#
function rev_length_scale(widths, ϕ_by_width, ϕ_bulk; tol = 1.0)
    idx = findlast(i -> abs(ϕ_by_width[i] - ϕ_bulk) > tol, eachindex(widths))
    idx === nothing && return widths[1]
    return widths[min(idx + 1, length(widths))]
end

rev_tol = 0.2 # tolerance percentage points around the bulk ϕ
rev_top_px = rev_length_scale(widths, ϕ_top_by_width, ϕ_top; tol = rev_tol)
rev_bot_px = rev_length_scale(widths, ϕ_bot_by_width, ϕ_bot; tol = rev_tol)

@info(

    "Estimated REV length scale (window size beyond which ϕ stays within tolerance $(rev_tol)pp of bulk ϕ):"
    ,
    REV_top_px = rev_top_px,
    REV_bot_px = rev_bot_px,
)

vlines!(
    ax_rev,
    [rev_top_px],
    color = top_color,
    linewidth = 3,
    linestyle = :dash,
    label = "REV Top $(rev_top_px) px",
)
vlines!(
    ax_rev,
    [rev_bot_px],
    color = bot_color,
    linewidth = 3,
    linestyle = :dash,
    label = "REV Bottom $(rev_bot_px) px",
)

delete!(rev_legend)
axislegend(ax_rev, position = :cb)
save("figures/berea_sandstone_rev.png", rev)
display(rev) # view the plot with REV markers
