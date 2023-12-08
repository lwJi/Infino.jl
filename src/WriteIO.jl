module WriteIO

using Plots

using DelimitedFiles

function dump(gfs, it)
  fname = "data/u.it" * lpad(string(it), 6, '0') * ".tsv"

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
end

function generate_gifs(dir_path)
  println("Generating gifs...")
  if isdir(dir_path)
    files = readdir(dir_path)
    a_psi = Animation()
    #a_Pi  = Animation()

    for fname in files
      file_path = joinpath(dir_path, fname)
      if isfile(file_path)
        println("Reading file: $file_path")
        data = readdlm(file_path, Float64, comments=true)
        # println(data[:, 5])

        x = []
        psi = []
        Pi = []
        for r in eachrow(data)
          if r[3] == 1
            push!(x, r[5])
            push!(psi, r[6])
            push!(Pi, r[7])
          end
        end
        plt_psi = plot(x, psi, ylim=(-1,1), label="psi")

        levs = Int(maximum(data[:, 3]))
        for l in 1:levs
          x = []
          psi = []
          Pi = []
          for r in eachrow(data)
            if r[3] == l
              push!(x, r[5])
              push!(psi, r[6])
              push!(Pi, r[7])
            end
          end
          plt_psi = scatter!(x, psi)
        end

        #plt_Pi  = plot(x, Pi, ylim=(-4,4), label="psi")
        #  plt_Pi  = scatter!(x, Pi)

        frame(a_psi, plt_psi)
        #frame(a_Pi, plt_Pi)

      end
    end

    gif(a_psi, "psi.gif")
    #gif(a_Pi, "Pi.gif")
  else
    println("Error: directory '$dir_path' does not exist!")
  end
end

end
