using DelimitedFiles

function read_bn_dat(file_path::String; width = 400, height = 400, T = UInt8)
    raw = readdlm(file_path)
    @assert size(raw) == (height, width) "expected $(height)x$(width), got $(size(raw))"
    matrix = T.(raw)
    void_space = findall(x -> x == false, matrix)
    solid_space = findall(x -> x == true, matrix)
    void_percentage = length(void_space) / (width * height) * 100
    @info("$file_path:",
        void = length(void_space),
        void_percentage = "$(round(void_percentage, digits=2))%",
        solid = length(solid_space),
        solid_percentage = "$(round(100 - void_percentage, digits=2))%"
    )

    return matrix
end
