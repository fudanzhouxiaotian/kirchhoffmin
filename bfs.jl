using LinearAlgebra
using Laplacians
using SparseArrays
using Arpack

function bfs(dmax,G)
    ###BFS
    n=G.n;
    f=1;
    r=0;
    b=zeros(1,n+1);
    b[1]=dmax;
    vis=zeros(1,n);
    vis[dmax]=1;
    len=zeros(1,n);
    len[dmax]=0;
    lenindex=[];

    while r<f
        r+=1;
        tt=Int(b[r]);
        nbr1=G.nbr[tt];
        for i=1:length(nbr1)
            t=nbr1[i];
            if vis[t]==0
                f+=1;
                b[f]=t;
                vis[t]=1;
                len[t]=len[tt]+1;
                push!(lenindex,t);
            end
        end
    end

    return len,lenindex;
end
