using Plots
matrix = readdlm("sparse.out")
heatmap(matrix, aspect_ratio=1.0)
savefig("picture.png")
