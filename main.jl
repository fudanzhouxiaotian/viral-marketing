include("graph.jl")
include("core.jl")
include("samp.jl")

using LinearAlgebra
using Graphs
# using Laplacians

fname = open("filename.txt", "r")
str   = readline(fname);
nn    = parse(Int, str);

type     = 1;  #### 1: directed graphs; 2: undirected graphs->directed graphs; 3:synthetic graphs
tot_test = 20; #### repeat times
k_max    = 50; #### maximum selected red nodes

for nnnn = 1:nn
    if type==1
        str = readline(fname);
        G   = get_graph_direct(str[1])
    elseif type==2
        str = readline(fname);
        G   = get_graph_undirect(str[1])
    elseif type==3
        # G   = BA(,,)  ### synthetic graphs will be written later
    end
    GG = turn_type(G);
    n = G.n;
    m = G.m;
    b_num_set = [2;5;10;20]
    pagerank_score = pagerank(GG);
    between_score  = betweenness_centrality(GG);
    close_score    = closeness_centrality(GG);
    indegree_score = indegree_centrality(GG);
    outdegree_score = outdegree_centrality(GG);
    for b_num in b_num_set
        #  random select blue nodes
        b_node_set = rand(1:n,2*b_num);
        b_node_set = union(b_node_set)[1:b_num];

        # baselines select red nodes
        tmp_pagerank_score = zeros(n);
        tmp_between_score  = zeros(n);
        tmp_close_score    = zeros(n);
        tmp_indegree_score = zeros(n);
        tmp_outdegree_score= zeros(n);

        
        tmp_pagerank_score .= pagerank_score;
        tmp_between_score  .= between_score;
        tmp_close_score    .= close_score;
        tmp_indegree_score .= indegree_score;
        tmp_outdegree_score.= outdegree_score;

        for b_select_node in b_node_set
            tmp_pagerank_score[b_select_node] = 0;
            tmp_between_score[b_select_node]  = 0;
            tmp_close_score[b_select_node]    = 0;
            tmp_indegree_score[b_select_node] = 0;
            tmp_outdegree_score[b_select_node]= 0;
        end
        r_node_set_1 = zeros(k_max);
        r_node_set_2 = zeros(k_max);
        r_node_set_3 = zeros(k_max);
        r_node_set_4 = zeros(k_max);
        r_node_set_5 = zeros(k_max);
        for k = 1:k_max
            push!(r_node_set_1, argmax(tmp_pagerank_score))
            tmp_pagerank_score[argmax(tmp_pagerank_score)]=0;
            push!(r_node_set_2, argmax(tmp_between_score))
            tmp_pagerank_score[argmax(tmp_between_score)]=0;
            push!(r_node_set_3, argmax(tmp_close_score))
            tmp_pagerank_score[argmax(tmp_close_score)]=0;
            push!(r_node_set_4, argmax(tmp_indegree_score))
            tmp_pagerank_score[argmax(tmp_indegree_score)]=0;
            push!(r_node_set_5, argmax(tmp_outdegree_score))
            tmp_pagerank_score[argmax(tmp_outdegree_score)]=0;
        end

        tot_time_1   = zeros(k_max);
        tot_red_1    = zeros(k_max);
        tot_time_2   = zeros(k_max);
        tot_red_2    = zeros(k_max);
        tot_time_3   = zeros(k_max);
        tot_red_3    = zeros(k_max);
        tot_time_4   = zeros(k_max);
        tot_red_4    = zeros(k_max);
        tot_time_5   = zeros(k_max);
        tot_red_5    = zeros(k_max);


        ### greedy algorithm to be written
        for testtimes = 1:tot_test    
            for k = 1:k_max
                t,r = final_state(G,b_node_set,r_node_set_1[1:k])
                tot_time_1[k] +=t/tot_test;
                tot_red_1[k]  +=r/tot_test;
                t,r = final_state(G,b_node_set,r_node_set_2[1:k])
                tot_time_2[k] +=t/tot_test;
                tot_red_2[k]  +=r/tot_test;
                t,r = final_state(G,b_node_set,r_node_set_3[1:k])
                tot_time_3[k] +=t/tot_test;
                tot_red_3[k]  +=r/tot_test;
                t,r = final_state(G,b_node_set,r_node_set_4[1:k])
                tot_time_4[k] +=t/tot_test;
                tot_red_4[k]  +=r/tot_test;
                t,r = final_state(G,b_node_set,r_node_set_5[1:k])
                tot_time_5[k] +=t/tot_test;
                tot_red_5[k]  +=r/tot_test;
            end
        end

        fout = open("ans.txt","a")
        println(fout,"b_num=",b_num,", baselines=page,bet,clo,indeg,outdeg,....")
        for k = 1:k_max
            println(fout,k,' ',tot_time_1[k],' ',tot_time_2[k],' ',tot_time_3[k],' ',tot_time_4[k],' ',tot_time_5[k],' ',tot_time_6[k])
        end
        for k = 1:k_max
            println(fout,k,' ',tot_red_1[k],' ',tot_red_2[k],' ',tot_red_3[k],' ',tot_red_4[k],' ',tot_red_5[k],' ',tot_red_6[k])
        end
        close(fout)
        
            
end
close(fname)
