struct Graph_direct
    n :: Int # |V|
    m :: Int # |E|
    u :: Array{Int32, 1}
    v :: Array{Int32, 1} # uv is an edge
#    w :: Array{Float64, 1} # weight of each edge
    nbr_in  :: Array{Array{Int32, 1}, 1}
    nbr_out :: Array{Array{Int32, 1}, 1}
end

include("core.jl")
using Graphs


function get_graph_direct(ffname)
    n = 0
    Label = Dict{Int32, Int32}()
    Origin = Dict{Int32, Int32}()
    #E = Set{Tuple{Int32, Int32, Float32}}()

    getID(x :: Int) = haskey(Label, x) ? Label[x] : Label[x] = n += 1

    fname = string("data/",ffname)
    fin = open(fname, "r")


    str = readline(fin)
    str = split(str)
    #n   = parse(Int, str[1])
    m   = parse(Int, str[3])
    u = Int32[]
    v = Int32[]

    tot = 0
    for i = 1 : m
        str = readline(fin)
        str = split(str)
        x   = parse(Int, str[1])
        y   = parse(Int, str[2])
        if x!=y
            u1 = getID(x)
            v1 = getID(y)
            Origin[u1] = x
            Origin[v1] = y
            push!(u, u1)
            push!(v, v1)
            tot += 1
        end
    end
    nbr_in=[ Int32[ ] for i in 1:n ]
    nbr_out=[ Int32[ ] for i in 1:n ]
    for i=1:tot
        u1=u[i];
        v1=v[i];
        push!(nbr_out[u1],v1);
        push!(nbr_in[v1],u1);
    end

    close(fin)
    return Graph_direct(n, tot, u, v,nbr_in,nbr_out)
end



function get_graph_undirect(ffname)
    n = 0
    Label = Dict{Int32, Int32}()
    Origin = Dict{Int32, Int32}()
    #E = Set{Tuple{Int32, Int32, Float32}}()

    getID(x :: Int) = haskey(Label, x) ? Label[x] : Label[x] = n += 1

    fname = string("data/",ffname)
    fin = open(fname, "r")


    str = readline(fin)
    str = split(str)
    #n   = parse(Int, str[1])
    m   = parse(Int, str[3])
    u = Int32[]
    v = Int32[]

    tot = 0
    for i = 1 : m
        str = readline(fin)
        str = split(str)
        x   = parse(Int, str[1])
        y   = parse(Int, str[2])
        if x!=y
            u1 = getID(x)
            v1 = getID(y)
            Origin[u1] = x
            Origin[v1] = y
            push!(u, u1)
            push!(v, v1)
            push!(u, v1)
            push!(v, u1)
            tot += 2
        end
    end
    nbr_in=[ Int32[ ] for i in 1:n ]
    nbr_out=[ Int32[ ] for i in 1:n ]
    for i=1:tot
        u1=u[i];
        v1=v[i];
        push!(nbr_out[u1],v1);
        push!(nbr_in[v1],u1);
    end

    close(fin)
    return Graph_direct(n, tot, u, v,nbr_in,nbr_out)
end


function turntype(G)
    g = SimpleDiGraph(G.n)
    for i = 1:G.m
        add_edge!(g,G.u[i],G.v[i])
    end
    return g;
end