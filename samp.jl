using SparseArrays
using Graphs
include("graph.jl")

function final_state(G,b_node_set,r_node_set)
    n = G.n;
    m = G.m;
    r = 0;
    round = 0;
    # uncolor = union(1:n)
    color   = zeros(n);
    color[b_node_set] .= -1;
    color[r_node_set] .= 1;
    # setdiff!(uncolor, b_node_set, r_node_set);
    neighb_is_col = Int32[];
    for i in b_node_set
        for j in G.nbr_in[i]
            if (color[j] == 0)
                push!(neighb_is_col, j);
            end
        end
    end
    union!(neighb_is_col);
    while size(neighb_is_col)[1]!=0
        round  += 1;
        tmp_col = spzeros(n);
        for i in neighb_is_col
            tmp_col[i] = color[rand(G.nbr_out[i])]
            if tmp_col[i]==1
                r += 1;
            end
        end
        color += tmp_col;
        
        # tmp_neighb = Int32[];
        for i in neighb_is_col
            if color[i]!=0
                setdiff!(neighb_is_col,i)
                for j in G.nbr_in[i]
                    if color[j]==0
                        union!(neighb_is_col,j)
                    end
                end
            end
        end
    end
    return round, r;
end

