include("graph.jl")
include("core.jl")
include("samp.jl")

using LinearAlgebra
using Graphs
# using Laplacians

fname = open("filename.txt", "r")
str   = readline(fname);
nn    = parse(Int, str);

type     =1;  #### 1: directed graphs; 2: undirected graphs->directed graphs; 3:synthetic graphs
tot_test = 50; #### repeat times

ba = [100,200,500,1000,2000,5000]

for nnnn = 1:nn
    t0 = time()
    if type==1
        str = readline(fname);
        G   = get_graph_direct(str)
    elseif type==2
        str = readline(fname);
        G   = get_graph_undirect(str)
    elseif type==3
        # G   = BA(,,)  ### synthetic graphs will be written later 6 networks
        GG = barabasi_albert(ba[nnnn],3)
        # println(GG)
        G  = backtype(GG)
        # println(G)
        str = "ba"
    end
    # println(GG)
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


    n = G.n;
    m = G.m;
    # println(G.n,' ',G.m)
    fout = open("time.txt","a")
    println(fout)
    println(fout,str,' ',G.n,' ',G.m," graphtype: ",type," epeated times: ",tot_test)
    close(fout)
    tot_time = zeros(n);
    for i=1:n
        for testtimes = 1:tot_test
            t=0;r=0;tc=0;
            t,r = final_state(G,union(i),[])
            tot_time[i] +=t/tot_test;
        end
    end
      
    fout = open("time.txt","a")
    println(fout,"max= ",maximum(tot_time)," min= ",minimum(tot_time)," average= ",sum(tot_time)/n)
    println(fout)
    for i=1:n
        # println(fout,tot_time[i])
    end
    close(fout)
    
    
        
            
end
close(fname)
