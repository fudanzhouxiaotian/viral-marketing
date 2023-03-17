using LinearAlgebra
using SparseArrays
using Random
using Arpack

include("samp.jl")

function qpl(G,k,b_node_set);
	s=zeros(k);
	for i=1:k
		s[i]=k;
	end
	t,r=final_state(G,b_node_set,s);
	indx=k;
	while s[1]<=n-k+1
		tmp_t,tmp_r=getans(L,n,s,k);
		if tmp_r>r
			r=tmp_r;
		end
		s[k]+=1;
		if s[k]>n
			s[k]-=1;
			while (indx>=1) && (s[indx]==n+indx-k)
				indx-=1;
				if indx==0
					indx=1;
					break;
				end
			end
			s[indx]+=1;
			for i=indx+1:k
				s[i]=s[i-1]+1;
			end
			indx=k;
		end
	end
	return tmp;
end




function lap_direct(G)
    F = zeros(G.n, G.n);
    for i = 1 : G.m
        F[G.u[i], G.v[i]] -= 1
        #F[G.v[i], G.u[i]] -= 1
        F[G.u[i], G.u[i]] += 1
        #F[G.v[i], G.v[i]] += 1
    end
    return F
end



function lapsp_direct(G)
	d=zeros(G.n);
	for i=1:G.m
		x=G.u[i];y=G.v[i];
		d[x]+=1;
	end
	uu=zeros(G.n+G.m)
	vv=zeros(G.n+G.m)
	ww=zeros(G.n+G.m)
	a=[i for i=1:G.n]
	uu[1:G.m]=G.u;
	uu[G.m+1:G.m+G.n]=a;
	vv[1:G.m]=G.v;
	vv[G.m+1:G.m+G.n]=a;
	ww[1:G.m].=-1;
	ww[G.m+1:G.n+G.m]=d;
	return sparse(uu,vv,ww)
end


function adjsp_direct(G)
	return sparse(G.u,G.v,ones(G.m),G.n,G.n)
end
