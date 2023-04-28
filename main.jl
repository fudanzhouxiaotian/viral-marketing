include("graph.jl")
include("core.jl")
include("samp.jl")
using LinearAlgebra
using Graphs
fname = open("filename.txt", "r")
str   = readline(fname);
nn    = parse(Int, str);
type     =1;  #### 1: directed graphs; 2: undirected graphs->directed graphs; 3:synthetic graphs
tot_test = 300; #### repeat times
k_max    = 10; #### maximum selected red nodes
# ba = [200,500,1000,2000,5000]
ba = [300]

for nnnn = 1:nn
    t0 = time()
    if type==1
        str = readline(fname);
        G   = get_graph_direct(str)
    elseif type==2
        str = readline(fname);
        G   = get_graph_undirect(str)
    elseif type==3
        # G   = BA(,,)  ### synthetic graphs will be written later
        GG = barabasi_albert(ba[nnnn],3)
        G  = backtype(GG)
        str = "ba"
    end
    t1 = time()
    GG = turntype(G);
    tmp_GG = strongly_connected_components(GG)
    find_max = zeros(size(tmp_GG)[1])
    size_max = 0;
    max_comp = Int[];
    for i = 1:size(tmp_GG)[1]
        if size(tmp_GG[i])[1]>size_max
            size_max = size(tmp_GG[i])[1]
            max_comp=tmp_GG[i]
        end
    end
    oon = G.n;
    oom = G.m;
    G=get_new(G,max_comp)
    GG = turntype(G);

    t2 = time()
    n = G.n;
    m = G.m;
    # println(G.n,' ',G.m)
    fout = open("newans.txt","a")
    println(fout)
    println(fout,str,' ',oon,' ',oom,' ',G.n,' ',G.m," graphtype: ",type," epeated times: ",tot_test)
    close(fout)
    b_num_set = [10;20;50;100]
    # b_num_set = [10]

    pagerank_score = pagerank(GG);
    between_score  = betweenness_centrality(GG);
    close_score    = closeness_centrality(GG);
    indegree_score = indegree_centrality(GG);
    outdegree_score = outdegree_centrality(GG);

    lab, no_use = label_propagation(GG);
    # println(maximum(lab))



    t3 = time()
    fout = open("newans.txt","a")
    println(fout)
    println(fout,"getingraph ",t1-t0,' ',t2-t1,' ',t3-t2)
    close(fout)
    for b_num in b_num_set
        t4 = time()
        #  random select blue nodes
        b_node_set = rand(1:n,2*b_num);
        b_node_set = union(b_node_set)[1:b_num];
        # baselines select red nodes
        xx = maximum(lab)
        uncolor_node = zeros(xx)
        for i=1:G.n
            if !(i in b_node_set)
                uncolor_node[lab[i]]+=1
            end
        end
        r_node_set_commu = zeros(Int,k_max)
        num = 0;
        while num<k_max
            xx = argmax(uncolor_node)
            for i=1:G.n
                if (lab[i]==xx) && !(i in b_node_set) && !(i in r_node_set_commu)
                    num+=1
                    uncolor_node[xx]-=1;
                    r_node_set_commu[num] = i
                    break
                end
            end
        end

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
        r_node_set_1 = zeros(Int,k_max);
        r_node_set_2 = zeros(Int,k_max);
        r_node_set_3 = zeros(Int,k_max);
        r_node_set_4 = zeros(Int,k_max);
        r_node_set_5 = zeros(Int,k_max);
        for k = 1:k_max
            r_node_set_1[k]= argmax(tmp_pagerank_score);
            tmp_pagerank_score[argmax(tmp_pagerank_score)]=0;
            r_node_set_2[k]= argmax(tmp_between_score);
            tmp_between_score[argmax(tmp_between_score)]=0;
            r_node_set_3[k]= argmax(tmp_close_score);
            tmp_close_score[argmax(tmp_close_score)]=0;
            r_node_set_4[k]= argmax(tmp_indegree_score);
            tmp_indegree_score[argmax(tmp_indegree_score)]=0;
            r_node_set_5[k]= argmax(tmp_outdegree_score);
            tmp_outdegree_score[argmax(tmp_outdegree_score)]=0;
        end
        t5 = time()
        greedy_node  = zeros(Int, k_max)
        tot_red_gre  = zeros(k_max);
        uncolor_set = union(1:n)
        setdiff!(uncolor_set,b_node_set)
        for k = 1:k_max
            red_num    = zeros(n);
            for to_selc_node in uncolor_set
                for testtimes = 1:tot_test
                    t=0;r=0;tc=0;
                    t,r=  final_state(G,b_node_set,union(greedy_node[1:k-1],to_selc_node))
                    red_num[to_selc_node] += r/ tot_test
                end
            end
            k_select_node = argmax(red_num)
            tot_red_gre[k] = maximum(red_num)
            setdiff!(uncolor_set, k_select_node)
            greedy_node[k] = k_select_node
        end
        t6 = time()
        fast_greedy_node  = zeros(Int, k_max)
        tot_red_fast_gre  = zeros(k_max);
        tmp_outdegree_score.= outdegree_score;
        tmp_outdegree_score[b_node_set].=0;
        uncolor_set = Int[];
        for i=1:100
            union!(uncolor_set,argmax(tmp_outdegree_score))
            tmp_outdegree_score[argmax(tmp_outdegree_score)]=0;
        end
        for k = 1:k_max
            red_num    = zeros(n);
            for i= 1: Int(round(G.n/100)+1) 
                to_selc_node = uncolor_set[i]
                for testtimes = 1:tot_test
                    t=0;r=0;tc=0;
                    t,r=  final_state(G,b_node_set,union(greedy_node[1:k-1],to_selc_node))
                    red_num[to_selc_node] += r/ tot_test
                end
            end
            k_select_node = argmax(red_num)
            tot_red_fast_gre[k] = maximum(red_num)
            setdiff!(uncolor_set, k_select_node)
            fast_greedy_node[k] = k_select_node
        end
        t6_fast = time()


        tot_red_1    = zeros(k_max);
        tot_red_2    = zeros(k_max);
        tot_red_3    = zeros(k_max);
        tot_red_4    = zeros(k_max);
        tot_red_5    = zeros(k_max);
        tot_red_6    = zeros(k_max);
        tot_red_7    = zeros(k_max);
        tot_red_8    = zeros(k_max);
        min_red_1    = ones(k_max).*G.n;
        min_red_2    = ones(k_max).*G.n;
        min_red_3    = ones(k_max).*G.n;
        min_red_4    = ones(k_max).*G.n;
        min_red_5    = ones(k_max).*G.n;
        min_red_6    = ones(k_max).*G.n;
        min_red_7    = ones(k_max).*G.n;
        min_red_8    = zeros(k_max).*G.n;
        sd_1 = zeros(k_max);
        sd_2 = zeros(k_max);
        sd_3 = zeros(k_max);
        sd_4 = zeros(k_max);
        sd_5 = zeros(k_max);
        sd_6 = zeros(k_max);
        sd_7 = zeros(k_max);
        sd_8 = zeros(k_max);
        t7 = time()
        t=0;r=0;tc=0;
            
        for k = 1:k_max
            test_tmp_1 = zeros(tot_test)
            test_tmp_2 = zeros(tot_test)
            test_tmp_3 = zeros(tot_test)
            test_tmp_4 = zeros(tot_test)
            test_tmp_5 = zeros(tot_test)
            test_tmp_6 = zeros(tot_test)
            test_tmp_7 = zeros(tot_test)
            test_tmp_8 = zeros(tot_test)
            for testtimes = 1:tot_test
                t=0;r=0;tc=0;
                t,r = final_state(G,b_node_set,r_node_set_1[1:k])
                tot_red_1[k]  +=r/tot_test;
                if r<min_red_1[k]
                    min_red_1[k]=r
                end
                test_tmp_1[testtimes] = r;
                t=0;r=0;tc=0;
                t,r = final_state(G,b_node_set,r_node_set_2[1:k])
                tot_red_2[k]  +=r/tot_test;
                if r<min_red_2[k]
                    min_red_2[k]=r
                end
                test_tmp_2[testtimes] = r;
                t=0;r=0;tc=0;
                t,r = final_state(G,b_node_set,r_node_set_3[1:k])
                tot_red_3[k]  +=r/tot_test;
                if r<min_red_3[k]
                    min_red_3[k]=r
                end
                test_tmp_3[testtimes] = r;
                t=0;r=0;tc=0;
                t,r = final_state(G,b_node_set,r_node_set_4[1:k])
                tot_red_4[k]  +=r/tot_test;
                if r<min_red_4[k]
                    min_red_4[k]=r
                end
                test_tmp_4[testtimes] = r;
                t=0;r=0;tc=0;
                t,r = final_state(G,b_node_set,r_node_set_5[1:k])
                tot_red_5[k]  +=r/tot_test;
                if r<min_red_5[k]
                    min_red_5[k]=r
                end
                test_tmp_5[testtimes] = r;
                t=0;r=0;tc=0;
                t,r = final_state(G,b_node_set,greedy_node[1:k])
                tot_red_6[k]  +=r/tot_test;
                test_tmp_6[testtimes] = r;
                t=0;r=0;tc=0;
                t,r = final_state(G,b_node_set,fast_greedy_node[1:k])
                tot_red_7[k]  +=r/tot_test;
                test_tmp_7[testtimes] = r;
                t=0;r=0;tc=0;
                t,r = final_state(G,b_node_set,r_node_set_commu[1:k])
                tot_red_8[k]  +=r/tot_test;
                test_tmp_8[testtimes] = r;
                if r>min_red_8[k]
                    min_red_8[k]=r
                end
            end
            sd_1[k] = sqrt(sum((test_tmp_1.-tot_red_1[k]).^2)/(tot_test))
            sd_2[k] = sqrt(sum((test_tmp_2.-tot_red_2[k]).^2)/(tot_test))
            sd_3[k] = sqrt(sum((test_tmp_3.-tot_red_3[k]).^2)/(tot_test))
            sd_4[k] = sqrt(sum((test_tmp_4.-tot_red_4[k]).^2)/(tot_test))
            sd_5[k] = sqrt(sum((test_tmp_5.-tot_red_5[k]).^2)/(tot_test))
            sd_6[k] = sqrt(sum((test_tmp_6.-tot_red_6[k]).^2)/(tot_test))
            sd_7[k] = sqrt(sum((test_tmp_7.-tot_red_7[k]).^2)/(tot_test))
            sd_8[k] = sqrt(sum((test_tmp_8.-tot_red_8[k]).^2)/(tot_test))
        end
        t8 = time()
        fout = open("newans.txt","a")
        println(fout,"b_num=",b_num,", baselines=page,bet,clo,indeg,outdeg,greedy,fastgreedy")
        println(fout,"time used: ",t5-t4,' ',t6-t5,' ',t7-t6,' ',t8-t7)
        println(fout,"fastgreedy time:",t6_fast-t6," greedy time:",t6-t5)
        println(fout,"red number:")
        # for k = 1:k_max
        #     println(fout,k,' ',tot_red_1[k]/G.n,' ',tot_red_2[k]/G.n,' ',tot_red_3[k]/G.n,' ',tot_red_4[k]/G.n,' ',tot_red_5[k]/G.n,' ',tot_red_6[k]/G.n,' ',tot_red_7[k]/G.n)
        # end
        for k = 1:k_max
            println(fout,k,' ',min_red_1[k]/G.n,' ',min_red_2[k]/G.n,' ',min_red_3[k]/G.n,' ',min_red_4[k]/G.n,' ',min_red_5[k]/G.n,' ',tot_red_6[k]/G.n,' ',tot_red_7[k]/G.n,' ',min_red_8[k]/G.n)
        end
        println(fout,"standard deviation")
        for k = 1:k_max
            println(fout,k,' ',sd_1[k]/G.n,' ',sd_2[k]/G.n,' ',sd_3[k]/G.n,' ',sd_4[k]/G.n,' ',sd_5[k]/G.n,' ',sd_6[k]/G.n,' ',sd_7[k]/G.n,' ',sd_8[k]/G.n)
        end
        println(fout)
        close(fout)
    end   
end
close(fname)
