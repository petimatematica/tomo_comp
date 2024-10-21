using Plots

include("procedures.jl")

# A, m, x, k = Landweber("Data82.mat")

n = 5000

A, m, x, k, V = SteepestDescent("Data328.mat", maxk=n)
img = heatmap(x, color = :greys, legend=false)

# X = Array{Plots.Plot{Plots.GRBackend}}(undef, n+1)

# for i in 1:n+1
#     X[i] = heatmap(V[i], color = :greys, legend=false, xlabel="iterada $(i-1)")
#     savefig(X[i], "D:/Usuario/Documents/Wellington/CNMAC2025/x-ray raw data/walnut/frame $(i-1)")
# end