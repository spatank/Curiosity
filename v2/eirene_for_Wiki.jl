
using Printf
using Random
using LinearAlgebra
using Eirene
using MAT
using Statistics
using Plots
using Glob

function constant_probability(n, p)
    G = zeros(n, n)
    for i = 1:n
        for j = 1:i-1
            r = rand(1)[1]
            if r < p
                G[i, j] = 1
                G[j, i] = 1
            end
        end
    end
    node_order = 1:n
    return G, node_order
end

function make_weighted_from_order(G, node_order)
    # adapted from code by ASB
    # original at https://github.com/BassettLab/Reorderability_scripts
    reordered_G = G[node_order, node_order] # often unnecessary
    n = length(node_order) # number of nodes
    val_mat = ones(n, n)
    for col = 1:n
        val_mat[1:col, col] = val_mat[1:col, col] * col
        val_mat[col, 1:col] = val_mat[col, 1:col] * col
    end
    weighted_G = reordered_G .* val_mat
    # replace 0 weighted edges with the largest edge weight possible
    # this is equivalent to assigning these edges the worst rank possible
    replace!(weighted_G, 0 => 2 * n)
    # weighted_G[findall(A -> A .== 0, weighted_G)] .= 2 * n would work too
    weighted_G[diagind(weighted_G)] .= 0 # set diagonal to 0
    return weighted_G # edges here are ranked by order of appearance
end

function bettiCurveFromBarcode(barcode_array,nNodes,nmats,maxDim)
    nNodes = Int(nNodes)
    nmats = Int(nmats)
    maxDim = Int(maxDim)
    bettiBar = zeros(nmats,maxDim)
    bettiCurve = zeros(nmats,nNodes+1,maxDim)
    birthCurve = zeros(nmats,nNodes,maxDim)
    deathCurve = zeros(nmats,nNodes,maxDim)
    for dimn in collect(1:maxDim)
        dimn = Int(dimn)
        for matn in collect(1:nmats)
            matn = Int(matn)
            bb = 0
            currentCurve = barcode_array[matn,:]
            currentCurveDim = currentCurve[dimn]
            for barn in collect(1:size(currentCurveDim,1))
                # Add to birth curve
                birthCurve[matn,Int(currentCurveDim[barn,1]),dimn] = birthCurve[matn,Int(currentCurveDim[barn,1]),dimn] .+1
                if currentCurveDim[barn,2]>nNodes
                    bettiCurve[matn,Int(currentCurveDim[barn,1]):Int(nNodes+1),dimn] = bettiCurve[matn,Int(currentCurveDim[barn,1]):Int(nNodes+1),dimn] .+1
                    bb = bb+(nNodes+1-currentCurveDim[barn,1])
                else
                    bettiCurve[matn,Int(currentCurveDim[barn,1]):Int(currentCurveDim[barn,2]),dimn] = bettiCurve[matn,Int(currentCurveDim[barn,1]):Int(currentCurveDim[barn,2]),dimn].+1
                    deathCurve[matn,Int(currentCurveDim[barn,2]),dimn] = deathCurve[matn,Int(currentCurveDim[barn,2]),dimn] .+1
                    bb = bb+(currentCurveDim[barn,2] - currentCurveDim[barn,1])
                end
            end
            bettiBar[matn,dimn] = deepcopy(bb)
        end
    end
    return bettiCurve, birthCurve, deathCurve, bettiBar
end

function plotBarcode(allPIs,nNodes,graphN,maxDim,fontSize)
    nNodes = Int(nNodes)
    graphn = Int(graphN)
    maxDim = Int(maxDim)
    counter1 = 0
    pbar = plot(1:6,zeros(6),c=:black)
    colors = [:blue :green :red]
    for dim in collect(1:maxDim)
        barn = barcode_array[graphN, dim]
        barn = barn[sortperm(barn[:,1]),:]
        nbars = size(barn)[1]
        for cntr1 in collect(1:nbars)
            birth = barn[cntr1,1]
            death = barn[cntr1,2]
            plot!([birth, death],[cntr1+counter1, cntr1+counter1],c=colors[dim], legend = false,
                            xlim = (0,nNodes), ytickfont = font(fontSize), xtickfont = font(fontSize))
        end
        display(pbar)
        counter1 = counter1+nbars
    end
    return pbar
end

pwd()

data_in_path = "/Users/sppatankar/Developer/My Passport/Curiosity/v2/Data/Wiki/Wiki_preprocessed_Eirene"
data_out_path = "/Users/sppatankar/Developer/My Passport/Curiosity/v2/Data/Wiki/Wiki_processed_Eirene"
all_topics = glob("*.mat", data_in_path)

topic = all_topics[6]

n_dims = 3 # number of dimensions to track persistent homology in

i = 6;
@printf("Topic: %s\n", all_topics[i])
all_vars = matread(all_topics[i])
nodes = all_vars["nodes"]
n = length(nodes)
adj = Array(all_vars["adj"])
edge_info = all_vars["edge_info"]
topic_ID = all_vars["topic"]
weighted_G = make_weighted_from_order(adj, 1:n)

save_string = string(data_out_path, "/", "weight_G_", topic_ID, ".mat")
matwrite(save_string,
    Dict("n" => n,
        "adj" => adj,
        "weighted_adj" => weighted_G))

barcode_array = Array{Array{Float64}}(undef, 3)
betti_curves = Array{Array{Float64}}(undef, 3)
C = Eirene.eirene(weighted_G, model = "vr", maxdim = n_dims, record = "none")
for k in 1:n_dims
    barcode_array[k] = barcode(C, dim = k)
    betti_curves[k] = betticurve(C, dim = k)
end
save_string = string(data_out_path, "/", "biophysics_Eirene.mat")
matwrite(save_string,
    Dict("n" => n,
        "adj" => adj,
        "weighted_adj" => weighted_G,
        "barcode_array" => barcode_array,
        "betti_curves" => betti_curves))
