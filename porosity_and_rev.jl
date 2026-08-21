#=
    Author: Viet M. Bui
    Date: 2026-08-21
    GLY 6826 - Hydrogeologic Modeling
    Homework 1: Porosity and Representative Elementary Volume (REV)

    Note: in order to fun the file, make sure to check for `data` 
    and `funcs` folders in the same dir
=#
include("funcs/read_bn_dat.jl")
using CairoMakie
using Colors

light_buff = colorant"#F0DC82"
buff = colorant"#70673b"

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
solid_volume_top = length(findall(x -> x == true, top_cross))
ϕ_top = void_volume_top / total_volume * 100
void_volume_bot = length(findall(x -> x == false, bot_cross))
solid_volume_bot = length(findall(x -> x == true, bot_cross))
ϕ_bot = void_volume_bot / total_volume * 100

@info(
    "Estimate ϕ of the rock for each cross section:",
    ϕ_Top = ϕ_top,
    ϕ_Bottom = ϕ_bot
)

#=
    c) For each cross section,

    calculate the average porosity within a square window centered on the image as a function of window size. Plot the results for both cross sections on the same graph.
=#
