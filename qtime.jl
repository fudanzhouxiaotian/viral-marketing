include("graph.jl")
include("core.jl")
include("samp.jl")
include("bfs.jl")

using LinearAlgebra
using Graphs
# using Laplacians

fname = open("filename.txt", "r")
str   = readline(fname);
nn    = parse(Int, str);

type     =2;  #### 1: directed graphs; 2: undirected graphs->directed graphs; 3:synthetic graphs
tot_test = 50; #### repeat times
qset= [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
ba = [300,500,1000,2000]
println(11)

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
    # println(GG)
    tmp_GG = strongly_connected_components(GG)
    # println(tmp_GG)
    find_max = zeros(size(tmp_GG)[1])
    size_max = 0;
    max_comp = Int[];
    for i = 1:size(tmp_GG)[1]
        if size(tmp_GG[i])[1]>size_max
            size_max = size(tmp_GG[i])[1]
            max_comp=tmp_GG[i]
        end
    end
    G=get_new(G,max_comp)
    GG = turntype(G);
    dd = zeros(G.n)
    # for i=1:G.n
    #     dd[i] = size(G.nbr_in[i])[1]
    # end
    # dmax = argmax(dd)
    # len,lenindex = bfs(dmax,G)
    # println(maximum(dd),' ',maximum(len))

    n = G.n;
    m = G.m;
    # println(G.n,' ',G.m)
    fout = open("qtime.txt","a")
    println(fout)
    println(fout,str,' ',G.n,' ',G.m," graphtype: ",type," repeated times: ",tot_test)
    println(str)
    close(fout)
    for q in qset
        # println(q)
        t_max = 0;
        t_min = 100000;
        t_tot = 0;
        for test = 1:tot_test
            b_node_set = Int[]
            for i = 1:n
                xx = rand()
                if xx < q
                    push!(b_node_set,i)
                end
            end
            # println(b_node_set)
            t_tmp ,r_tmp = final_state(G,b_node_set,Int[])
            if t_tmp>t_max
                t_max = t_tmp
            end
            if t_tmp<t_min
                t_min = t_tmp
            end
            t_tot+=t_tmp/tot_test
        end            
        fout = open("qtime.txt","a")
        # println(fout,"q=",q," max= ",t_max," min= ",minimum(tot_time)," average= ",sum(tot_time)/n)
        println(fout, q,' ',t_max,' ',t_min,' ',t_tot)
        # println( q,' ',t_max,' ',t_min,' ',t_tot)
        # println(fout)
        close(fout)
    end
    # close(fout)   
end
close(fname)
