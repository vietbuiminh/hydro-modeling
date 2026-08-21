function read_bn_dat(file_path::String; width = 400, height = 400, T = UInt8)
    matrix = Matrix{T}(undef, width, height)
    open(file_path, "r") do io
        read!(io, matrix)
    end
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
