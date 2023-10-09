using LinearAlgebra
n=3;
a = zeros(n,n)
for i=1:n
    for j=i+1:n
        a[i,j]=rand([1,-1,0])
        a[j,i]=a[i,j]
    end
end
for i=1:n
    for j=1:n
        if i!=j
            a[i,i]+=-a[i,j]
        end
    end
end


function pp(a,n)
    for i=1:n
        for j=1:n
            print((a[i,j]),' ')
        end
        println()
    end    
    println()
end

a = [0 -1 1;-1 0 1;1 1 -2]
pp(a,n)

J=ones(n,n)./n;

if rank(a)==n-1
    println("inv:")
    b = inv(a+J)-J 
    # pp(inv(a+J)-J,n)
    global c = b;
    for i=1:3
        global c = c*b;
        pp(b^i,n)
        # pp(b^(1/i),n)
    end
end

# global b=a*a
# pp(b,n)
# a=b;

# println("new")
# pp(inv(b+J)-J,n)

for i=1:3
    # global b=b*a
    pp(a^i,n)
    # pp(a^(1/i),n)
end