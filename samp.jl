using SparseArrays
using Graphs
include("graph.jl")

function final_state(G,b_node_set,r_node_set)
    n = G.n;
    m = G.m;
    rr = size(r_node_set)[1];
    round = 0;
    tot_color = size(b_node_set)[1]+size(r_node_set)[1];
    color   = zeros(n);
    color[b_node_set] .= -1;
    color[r_node_set] .= 1;
    while tot_color<n
        round+=1;
        for i=1:n
            if color[i]==0
                color[i]=color[rand(G.nbr_out[i])];
                if color[i]==1
                    rr+=1
                end
                if color[i]!=0
                    tot_color+=1
                end
            end
        end
    end
    return round, rr;
end



function final_state_complex(G,b_node_set,r_node_set)
    n = G.n;
    m = G.m;
    r = size(r_node_set)[1];
    round = 0;
    tot_color = size(b_node_set)[1]+size(r_node_set)[1];
    # println("**********")
    # println(b_node_set);
    # println(r_node_set)
    # uncolor = union(1:n)
    color   = zeros(n);
    color[b_node_set] .= -1;
    color[r_node_set] .= 1;
    # println(color)
    # setdiff!(uncolor, b_node_set, r_node_set);
    neighb_is_col = Int32[];
    for i in union(b_node_set,r_node_set)
        for j in G.nbr_in[i]
            if (color[j] == 0)
                union!(neighb_is_col, j);
            end
        end
    end
    # union!(neighb_is_col);
    # println(neighb_is_col)
    while size(neighb_is_col)[1]!=0
        round  += 1;
        tmp_col = zeros(n);
        new_col = Int[];
        for i in neighb_is_col
            tmp_col[i] = color[rand(G.nbr_out[i])]
            if tmp_col[i]!=0
                tot_color += 1;
                union!(new_col,G.nbr_in[i])
            end
            if tmp_col[i]==1
                r += 1;
            end
        end
        # println(tmp_col)
        # println(r)
        color += tmp_col;
        # println(color)
        
        # tmp_neighb = Int32[];
        not_color = Int32[];
        
        union!(neighb_is_col,new_col);
        for i in neighb_is_col
            if color[i]!=0
                union!(not_color,i)
                # setdiff!(neighb_is_col,i)
            end
        end
        setdiff!(neighb_is_col,not_color)
        # if size(neighb_is_col)[1]+tot_color!=16
            # println(size(neighb_is_col)[1]+tot_color)
        # end
    end
    return round, r, tot_color;
end

