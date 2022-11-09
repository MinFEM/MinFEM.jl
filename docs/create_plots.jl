using MinFEM, Plots

path = "src/assets/fem/"

fmt = "svg" # also "png", "pdf"

pyplot(leg=false,linewidth=1,ticks=true)
plot([0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], label="T",
        title="Reference Element", size = (350, 300))
savefig(path * "reference_element.$fmt")

pyplot(leg=false,linewidth=1,ticks=false)
p1 = surface([0.0, 1.0, 0.0], [0.0, 0.0, 1.0],[0,0,0],
                colorbar = :none, color=:viridis)
surface!([0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0,0.0],
            colorbar = :none, color=:viridis)
p2 = surface([0.0, 1.0, 0.0], [0.0, 0.0, 1.0],[0,0,0],
                colorbar = :none, color=:viridis)
surface!([0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0,0.0],
            colorbar = :none, color=:viridis)
p3 = surface([0.0, 1.0, 0.0], [0.0, 0.0, 1.0],[0,0,0],
                colorbar = :none, color=:viridis)
surface!([0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, 1.0],
            colorbar = :none, color=:viridis)
p4 = plot(p1, p2, p3, layout=(1,3),size=(900,350))
savefig(path * "local_basis.$fmt")

pyplot(leg=false,linewidth=1,ticks=false)
p1 = surface([0.0, -0.2, 1.0, 1.2, 0.4, 0.5], [0.0, 0.5, 0.0, 1.0, 1.2, 0.5],
                [0,0,0,0,0,1.0], camera=(80,30), colorbar=:none, color=:viridis)
p2 = surface([0.0, -0.2, 1.0, 1.2, 0.4, 0.5], [0.0, 0.5, 0.0, 1.0, 1.2, 0.5],
                [0,0,0,0,0,1.0], camera=(0,90), colorbar=:none, color=:viridis)
p3 = plot(p1, p2, layout=(1,2), size=(900,300))
savefig(path * "global_basis.$fmt")

quadX = quadrature_points(2, 3)

quadW = quadrature_weights(2, 3)

pyplot(leg=false,linewidth=1,ticks=true)
plot([0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], label="T", title="Quadrature Formula")
scatter!([x[1] for x in quadX], [x[2] for x in quadX], zcolor=quadW,
            label="Quadrature nodes", size = (350, 300),colorbar=false)
savefig(path * "quadrature.$fmt")

pyplot(leg=false,linewidth=1,ticks=true)
p1 = plot([0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], label="T", title="Reference Element")
p2 = plot([0.3, 0.8, 0.6, 0.3], [0.3, 0.4, 0.9, 0.3], label="T", title="Physical Element",
            xlims=(0,1), ylims=(0,1))
plot(p1, p2, layout = (1,2), size=(600,250))
savefig(path * "physical_element.$fmt")