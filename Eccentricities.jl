include("graph.jl")
include("core.jl")
include("bfs.jl")
using LinearAlgebra
using Laplacians
using SparseArrays
using Arpack


       
function Ecc_Comp(G)
    L=lapsp(G);
    n=G.n;
    m=G.m;
        degr=0;
        dmax=0;
        for i=1:n
            if L[i,i]>degr
                degr=L[i,i];
                dmax=i;
            end
        end
        ecc=zeros(1,n);#存ecc值
        ecc2=zeros(1,n);#存下标
        ecclow=zeros(1,n);
        eccup=zeros(1,n);
        len,lenindex=bfs(dmax,G);
        eccindex=argmax(len)[2];
        ecc2[dmax]=eccindex;
        ecc[dmax]=len[eccindex];
        for i=1:n
            if i!=dmax
                ecclow[i]=max(len[i],ecc[dmax]-len[i]);
                eccup[i]=ecc[dmax]+len[i];
            end
        end
        k1=length(lenindex);
        k2=100;
        for j=1:k2
            v1=lenindex[k1+1-j];
            lenv1,lenindexv1=bfs(v1,G);
            for i=1:n
                # if i!=dmax
                    ecclow[i]=max(ecclow[i],lenv1[i]);
                    eccup[i]=min(eccup[i],max(ecclow[i],len[v1]+len[i]));
                    if ecclow[i]==eccup[i]
                        ecc2[i]=v1;
                    end
                # end
            end
        end
    # end
    
    return ecclow;
end
