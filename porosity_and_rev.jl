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
│   void = 95018
│   void_percentage = "59.39%"
│   solid = 64982
└   solid_percentage = "40.61%"
┌ Info: data/berea_xsection_bot.dat:
│   void = 93357
│   void_percentage = "58.35%"
│   solid = 66643
└   solid_percentage = "41.65%"
=#

#= 
    a) Plot the two micro-CT cross sections.
=#
sandx = Figure(size = (500, 1000), dpi = 300, title = "Berea Sandstone XSection")
ax_top = Axis(sandx[1, 1], title = "Top")
ax_bot = Axis(sandx[2, 1], title = "Bottom")

binary = top_cross .>= true
heatmap!(ax_top, binary, colormap = [:white, :black])

Label(
    sandx[0, :], text = "Berea Sandstone XSection",
    font = :bold,
    fontsize = 23, tellwidth = false,
)

sandx # view the plot

#=
    b) Estimate the porosity of the rock from each cross section.
=#

#=
    c) For each cross section,
    calculate the average porosity within a square window centered on the image as a function of window size. Plot the results for both cross sections on the same graph.
=#
