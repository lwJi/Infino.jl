module WriteIO

using Plots
using DelimitedFiles

function dump(dir_path, gfs, it)
  if isdir(dir_path)
    fname = dir_path * "/u.it" * lpad(string(it), 6, '0') * ".tsv"
    open(fname, "w") do file
      println(file, "# 1:iteation 2:time 3:level 4:i 5:x 6:psi 7:Pi")
      for l in 1:length(gfs.levs)
        lev = gfs.grid.levs[l]
        levfs = gfs.levs[l]
        for i in 1:length(levfs.x)
          println(file,
                  it,
                  " ", lev.time,
                  " ", l,
                  " ", i,
                  " ", levfs.x[i],
                  " ", levfs.u[1][i],
                  " ", levfs.u[2][i])
        end
      end
    end
  else
    println("Error: directory '$dir_path' does not exist!")
  end
end

function generate_gifs(dir_path)
  println("Generating gifs...")
  if isdir(dir_path)
    files = readdir(dir_path)
    tsv_files = filter(f -> endswith(f, ".tsv"), files)
    a_psi = Animation()
    a_Pi  = Animation()
    for fname in tsv_files
      file_path = joinpath(dir_path, fname)
      if isfile(file_path)
        println("Reading file: $file_path")
        data = readdlm(file_path, Float64, comments=true)
        # plot psi
        x   = data[findall(x -> x == 1, data[:, 3]), 5]
        psi = data[findall(x -> x == 1, data[:, 3]), 6]
        plt_psi = plot(x, psi, ylim=(-1,1), label="psi")
        for l in 1:Int(maximum(data[:, 3]))
          x   = data[findall(x -> x == l, data[:, 3]), 5]
          psi = data[findall(x -> x == l, data[:, 3]), 6]
          plt_psi = scatter!(x, psi)
        end
        frame(a_psi, plt_psi)
        # plot Pi
        x  = data[findall(x -> x == 1, data[:, 3]), 5]
        Pi = data[findall(x -> x == 1, data[:, 3]), 7]
        plt_Pi = plot(x, Pi, ylim=(-4,4), label="Pi")
        #for l in 1:Int(maximum(data[:, 3]))
        #  x  = data[findall(x -> x == l, data[:, 3]), 5]
        #  Pi = data[findall(x -> x == l, data[:, 3]), 7]
        #  plt_Pi = scatter!(x, Pi)
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

end
