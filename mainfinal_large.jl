include("graph.jl")
include("core.jl")
# include("schur.jl")
include("edgecore.jl")
include("bb.jl")
include("bfs.jl")
include("Eccentricities.jl")

using LinearAlgebra
using Laplacians

fname = open("filename.txt", "r")
str   = readline(fname);
nn     = parse(Int, str);
tbn = 100; ###凸包点数
eps = 0.5;
kmax=50;
t = 40;
for nnnn=1:nn
    str = readline(fname);
    str = split(str);
    G   = get_graph(str[1]);
    on=G.n;om=G.m;
    Gc=findconnect(G);
    G=Gc;
    n=G.n;m=G.m;
    

    fans1=open("ans.txt","a");
    println(fans1,str[1],' ',on,' ',om,' ',n,' ',m)
    close(fans1)


    L=lapsp(G)
    A=adjsp(G)

    #####exact_begin simplegreedy
    # t11=time()
    # J = ones(n,n)./n;
    # invL = inv(L+J)-J;
    # ans1x = zeros(kmax)
    # ans1y = zeros(kmax)
    # for k=1:kmax
    #     e1=Int[];
    #     e2=Int[];
    #     aa=Float32[];
    #     for i=1:n
    #         for j=i+1:n
    #             if L[i,j]==0
    #                 push!(e1,i)
    #                 push!(e2,j)
    #                 # e=zeros(n)
    #                 # e[i]=1;e[j]=-1;
    #                 push!(aa,norm(invL[i,:]-invL[j,:],2)^2/(1+invL[i,i]+invL[j,j]-2*invL[i,j]))
    #             end
    #         end
    #     end
    #     selcedge = argmax(aa)
    #     ans1x[k]=e1[selcedge]
    #     ans1y[k]=e2[selcedge]
    #     # L[e1[selcedge],e2[selcedge]]=-1;
    #     # L[e2[selcedge],e1[selcedge]]=-1;
    #     e=zeros(n)
    #     e[e1[selcedge]]=1;
    #     e[e2[selcedge]]=-1;
    #     invL = invL - (invL*e)*(e'*invL)/(1+e'*invL*e)
    # end
    # t12=time()
    # println("exact complete")
    ###exact_end


    ######## solver+jl begin
    # L=lapsp(G)
    # A=adjsp(G)
    # t21=time()
    # ans2x = zeros(kmax)
    # ans2y = zeros(kmax)
    # uu = zeros(m)
    # vv = zeros(m)
    # uu .= G.u;
    # vv .= G.v;
    # B = sparse(1:m,uu,ones(m),m,n)-sparse(1:m,vv,ones(m),m,n)
    # for k=1:kmax
    #     # t = Int(round(log(n)/eps^2))
    #     f = approxchol_lap(A)
    #     coor = zeros(n,t)
    #     for i=1:t
    #         xx = rand(Int8[1,-1],n)
    #         yy = f(xx)
    #         coor[1:n,i]=yy;
    #     end
    #     coor2 = zeros(n,t)
    #     for i=1:t
    #         xx = B'*rand(Int8[1,-1],m+k-1);
    #         coor2[1:n,i]=f(xx);
    #     end

    #     ###find optimal
    #     tmp1=0;tmp2=0;tmp3=0;
    #     for i=1:n
    #         for j=i+1:n
    #             # e=zeros(n)
    #             # e[conv[i]]=1;
    #             # e[conv[j]]=-1;
    #             tmpp=norm(coor[i,:]-coor[j,:])^2/(1+norm(coor2[i,:]-coor2[j,:])^2)
    #             if  tmpp>tmp3
    #                 tmp3=tmpp
    #                 tmp1=i;
    #                 tmp2=j;
    #             end
    #         end
    #     end
    #     ans2x[k]=tmp1;
    #     ans2y[k]=tmp2;
    #     A[tmp1,tmp2]=1;
    #     A[tmp2,tmp1]=1;
    #     push!(uu,tmp1)
    #     push!(vv,tmp2)
    #     B = sparse(1:(m+k),uu,ones(m+k),m+k,n)-sparse(1:(m+k),vv,ones(m+k),m+k,n)
    # end
    # t22=time()
    # ######### solver+jl end
    # println("solver complete")

    # ######## solver+jl+kooji begin
    # L=lapsp(G)
    # A=adjsp(G)
    # t31=time()
    # ans3x = zeros(kmax)
    # ans3y = zeros(kmax)
    # uu = zeros(m)
    # vv = zeros(m)
    # uu .= G.u;
    # vv .= G.v;
    # B = sparse(1:m,uu,ones(m),m,n)-sparse(1:m,vv,ones(m),m,n)
    # coor2 = zeros(n,t)
    #     for i=1:t
    #         xx = B'*rand(Int8[1,-1],m);
    #         coor2[1:n,i]=f(xx);
    #     end
    #     lii=zeros(n)
    # for i=1:n
    #     lii[i]=norm(coor2[i,:])
    # end
    # bbb = sort(lii)
    # gd = zeros(n)
    # for i=1:n
    #     if lii[i]>=bbb[Int(round(n-n/sqrt(kmax)))]
    #         gd[i]=1
    #     end
    # end

    # for k=1:kmax
    #     # t = Int(round(log(n)/eps^2))
    #     f = approxchol_lap(A)
    #     coor = zeros(n,t)
    #     for i=1:t
    #         xx = rand(Int8[1,-1],n)
    #         yy = f(xx)
    #         coor[1:n,i]=yy;
    #     end
    #     coor2 = zeros(n,t)
    #     for i=1:t
    #         xx = B'*rand(Int8[1,-1],m+k-1);
    #         coor2[1:n,i]=f(xx);
    #     end

    #     ###find optimal
    #     tmp1=0;tmp2=0;tmp3=0;
    #     for i=1:n
    #         if gd[i]==1
    #         for j=i+1:n
    #             if gd[j]==1
    #             # e=zeros(n)
    #             # e[conv[i]]=1;
    #             # e[conv[j]]=-1;
    #             tmpp=tmpp=norm(coor[i,:]-coor[j,:])^2/(1+norm(coor2[i,:]-coor2[j,:])^2)
    #             if  tmpp>tmp3
    #                 tmp3=tmpp
    #                 tmp1=i;
    #                 tmp2=j;
    #             end
    #         end
    #         end
    #         end
    #     end
    #     ans3x[k]=tmp1;
    #     ans3y[k]=tmp2;
    #     A[tmp1,tmp2]=1;
    #     A[tmp2,tmp1]=1;
    #     push!(uu,tmp1)
    #     push!(vv,tmp2)
    #     B = sparse(1:(m+k),uu,ones(m+k),m+k,n)-sparse(1:(m+k),vv,ones(m+k),m+k,n)
    # end
    # t32=time()
    # ######### solver+jl+kooji end
    # println("kooji complete")



    #######approx_fast gradfast  begin fast alg with approx marginal
 
    L=lapsp(G)
    A=adjsp(G)
    t41=time()
    ans4x = zeros(kmax)
    ans4y = zeros(kmax)
    for k=1:kmax
        # t = Int(round(log(n)/eps^2))
        f = approxchol_lap(A)
        coor = zeros(n,t)
        for i=1:t
            xx = rand(Int8[1,-1],n)
            yy = f(xx)
            coor[1:n,i]=yy;
        end
        ###convex hull
        conv = bb(coor,tbn)
        # println(conv)
        ###find optimal
        num_conv = length(conv)
        tmp1=0;tmp2=0;tmp3=0;
        for i=1:num_conv
            for j=i+1:num_conv
                # e=zeros(n)
                # e[conv[i]]=1;
                # e[conv[j]]=-1;
                tmpp=norm(coor[conv[i],1:t]-coor[conv[j],1:t],2)
                if  tmpp>tmp3
                    tmp3=tmpp
                    tmp1=conv[i];
                    tmp2=conv[j];
                end
            end
        end
        ans4x[k]=tmp1;
        ans4y[k]=tmp2;
        A[tmp1,tmp2]=1;
        A[tmp2,tmp1]=1;
    end
    t42=time()
    #######approx_fast_end  fastgrad
    println("gradfast complete")

    ######gradfast+ fast

        L=lapsp(G)
        A=adjsp(G)
        t51=time()
        ans5x = zeros(kmax)
        ans5y = zeros(kmax)
        # t = Int(round(log(n)/eps^2))
        Q = rand(Int8[1,-1],n,t)
        f = approxchol_lap(A)
        coor = zeros(n,t)
        for i=1:t
            # xx = rand(Int8[1,-1],n)
            yy = f(Q[1:n,i])
            coor[1:n,i]=yy;
        end
        # ecc=zeros(n)
        ecc = Ecc_Comp(G)
        outnodes = Int[];
        for i=1:G.n
            p=1;
            for j in G.nbr[i]
                if ecc[j]>=ecc[i]
                    p=0
                    break
                end
            end
            if p==1
                push!(outnodes,i)
            end
        end
        # println(length(outnodes))
       #  maxecc = maximum(ecc)
       #  tj = zeros(Int(maxecc))
       #  for i=1:n
       #      tj[Int(ecc[i])]+=1;
       #  end
       #  xx = argmax(tj)
       #  xx = min(xx,maxecc-1)
       #  for i=1:G.n
       #      if ecc[i]>xx
       #          push!(outnodes,i)
       #      end
       #  end
       #  println(length(outnodes))
    
        for k=1:kmax
            f = approxchol_lap(A)
            ###convex hull
            conv = bb(coor[outnodes,:],min(tbn,length(outnodes)))
            # println(conv)
            ###find optimal
            num_conv = length(conv)
            tmp1=0;tmp2=0;tmp3=0;
            for i=1:num_conv
                for j=i+1:num_conv
                    tmpp=norm(coor[outnodes[conv[i]],1:t]-coor[outnodes[conv[j]],1:t],2)
                    if  tmpp>tmp3
                        tmp3=tmpp
                        tmp1=outnodes[conv[i]];
                        tmp2=outnodes[conv[j]];
                    end
                end
            end
            ans5x[k]=tmp1;
            ans5y[k]=tmp2;
    
            e = zeros(n)
            e[tmp1]=1;
            e[tmp2]=-1;
            tmpvec = f(e)
            fm = (e'*tmpvec)[1];
            ####update node coor
            pp = Q'*tmpvec
            for cn in outnodes
                # cn = tt;
                coor[cn,:]-=pp*tmpvec[cn]/fm;
            end
            A[tmp1,tmp2]=1;
            A[tmp2,tmp1]=1;
        end
        t52=time()
    println("gradfast+ complete")
    #######convex prune fast end



    #######one convex hull
    L = lapsp(G)
    A = adjsp(G)
    t61=time()
    ans6x = zeros(kmax)
    ans6y = zeros(kmax)
    # t = Int(round(log(n)/eps^2))
    Q = rand(Int8[1,-1],n,t)
    f = approxchol_lap(A)
    coor = zeros(n,t)
    for i=1:t
        # xx = rand(Int8[1,-1],n)
        yy = f(Q[1:n,i])
        coor[1:n,i]=yy;
    end
    conv = bb(coor,tbn*2)
    for k=1:kmax
        ###find optimal
        f = approxchol_lap(A)
        num_conv = length(conv)
        tmp1=0;tmp2=0;tmp3=0;
        for i=1:num_conv
            for j=i+1:num_conv
                e=zeros(n)
                e[conv[i]]=1;
                e[conv[j]]=-1;
                tmpp=norm(coor[conv[i],1:t]-coor[conv[j],1:t],2)
                if  tmpp>tmp3
                    tmp3=tmpp
                    tmp1=conv[i];
                    tmp2=conv[j];
                end
            end
        end
        ans6x[k]=tmp1;
        ans6y[k]=tmp2;
        e = zeros(n)
        e[tmp1]=1;
        e[tmp2]=-1;
        # setdiff!(conv,tmp1,tmp2)
        tmpvec = f(e)
        fm = (e'*tmpvec)[1];
        ####update node coor
        qq = Q'*tmpvec;
        for cn in conv
            # cn = conv[tt];
            coor[cn,:]-=qq*tmpvec[cn]/fm;
        end
        A[tmp1,tmp2]=1;
        A[tmp2,tmp1]=1;
    end
    t62=time()
    #######one convex hull fast end
    println("oneconvex complete")



    #######convex prune end
    
    ans7x=rand(1:n,kmax)
    ans7y=rand(1:n,kmax)
    A=adjsp(G)
    ans8x=zeros(kmax)
    ans8y=zeros(kmax)
    for k=1:kmax
        l,r = eigcal(A)
        ans8x[k]=argmax(r)[1]
        ans8y[k]=argmin(r)[1]
        A[Int(ans8x[k]),Int(ans8y[k])]+=1;
        A[Int(ans8y[k]),Int(ans8x[k])]+=1;
    end
    
    fout = open("ans.txt","a")
    println(fout,t42-t41,' ',t52-t51,' ',t62-t61,' ')
    A = adjsp(G)

    t_fin1=time()
    iniK = getkir(A,n)
    kir4=zeros(kmax)
    kir4[1]=iniK+kir(ans4x,ans4y,1,A,n)
    kir5=zeros(kmax)
    kir5[1]=iniK+kir(ans5x,ans5y,1,A,n)
    kir6=zeros(kmax)
    kir6[1]=iniK+kir(ans6x,ans6y,1,A,n)
    kir7=zeros(kmax)
    kir7[1]=iniK+kir(ans7x,ans7y,1,A,n)
    kir8=zeros(kmax)
    kir8[1]=iniK+kir(ans8x,ans8y,1,A,n)
    for k=2:kmax
        kir4[k]=kir4[k-1]+kir(ans4x,ans4y,k,A,n)
        kir5[k]=kir5[k-1]+kir(ans5x,ans5y,k,A,n)
        kir6[k]=kir6[k-1]+kir(ans6x,ans6y,k,A,n)
        kir7[k]=kir7[k-1]+kir(ans7x,ans7y,k,A,n)
        kir8[k]=kir8[k-1]+kir(ans8x,ans8y,k,A,n)
    end
    for k=1:kmax
        println(fout,k,' ',kir4[k],' ',kir5[k],' ',kir6[k],' ',kir7[k],' ',kir8[k])
    end
    t_fin2=time()
    println(fout,t_fin2-t_fin1)
    println(fout)
    close(fout)
end

close(fname)