using Plots
using DelimitedFiles

function generate_gifs()
    if length(ARGS) < 1
        println("Usage: julia generate_gifs.jl data_path")
        exit(1)
    end

    println("Generating gifs...")
    dir_path = ARGS[1]
    if isdir(dir_path)
        files = readdir(dir_path)
        tsv_files = filter(f -> endswith(f, ".tsv"), files)
        a_psi = Animation()
        a_Pi = Animation()
        for fname in tsv_files
            file_path = joinpath(dir_path, fname)
            if isfile(file_path)
                println("Reading file: $file_path")
                data = readdlm(file_path, Float64, comments = true)
                # plot psi
                x = data[findall(x -> x == 1, data[:, 3]), 5]
                psi = data[findall(x -> x == 1, data[:, 3]), 6]
                plt_psi = plot(x, psi, ylim = (-1, 1), label = "psi")
                for l = 1:Int(maximum(data[:, 3]))
                    x = data[findall(x -> x == l, data[:, 3]), 5]
                    psi = data[findall(x -> x == l, data[:, 3]), 6]
                    plt_psi = scatter!(x, psi, label = "")
                end
                frame(a_psi, plt_psi)
                # plot Pi
                x = data[findall(x -> x == 1, data[:, 3]), 5]
                Pi = data[findall(x -> x == 1, data[:, 3]), 7]
                plt_Pi = plot(x, Pi, ylim = (-4, 4), label = "Pi")
                #for l in 1:Int(maximum(data[:, 3]))
                #  x  = data[findall(x -> x == l, data[:, 3]), 5]
                #  Pi = data[findall(x -> x == l, data[:, 3]), 7]
                #  plt_Pi = scatter!(x, Pi, label="")
                #end
                frame(a_Pi, plt_Pi)
            end
        end
        gif(a_psi, dir_path * "/psi.gif")
        gif(a_Pi, dir_path * "/Pi.gif")
    else
        println("Error: directory '$dir_path' does not exist!")
    end
end

generate_gifs()
