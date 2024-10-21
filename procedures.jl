using MAT, Plots, LinearAlgebra

#= 

Function: loadfile
Description: retrieve the matrices from the files: Data82.mat, Data164.mat and Data328.mat
Arguments: filename - name of the file to be loaded (String)

=#
function loadfile(filename)
    data = matread("D:/Usuario/Documents/Wellington/CNMAC2025/x-ray raw data/walnut/"*filename)
    A = data["A"]
    m = data["m"]
    return A, m
end

#= 

Function: testnorm
Description: test if the norm of a matrix is lesser than or equal to a fixed positive number
Arguments: A - matrix
           x - fixed positive number
=#

function testnorm(A, x)
    
    v = [1,2,Inf] # vector with possible norm types

    for i in 1:3
        if ~(norm(A, v[i])>x)
            return v[i]
        end
    end

    return 0

end

#= 

Function: standardform
Description: piles columns of a given matrix to form a column vector
Arguments: A - matrix

=#

function standardform(A)

    m = size(A)[1]
    n = size(A)[2]

    V = zeros(m*n)

    for j in 1:n
        for i in 1:m
            V[i+(j-1)*m]=A[i,j]
        end
    end

    return V

end

#= 

Function: SquareReconstruction
Description: turn a vector (with the number of entries being a perfect square) into a square matrix
Arguments: V - vector

=#

function SquareReconstruction(V)

    n = Int64.(sqrt(size(V)[1]))
    m = zeros(n,n)

    for j in 1:n
        for i in 1:n 
            m[i,j]=V[i+(j-1)*n]
        end
    end

    return(m)

end

#= 

Function: Landweber
Description: implements the Landweber's iterative regularization method for inverse problems
Arguments: filename - name of the file to be loaded (String)
Opitional arguments: tol - tolerance given
                     maxk - maximum number of iterates allowed

=#

function Landweber(filename; tol=1.e-4, maxk=100)

    A, m = loadfile(filename)  # retrieving the matrices
    n = testnorm(A, 1) # because the convergence is done assuming that ||A|| is lesser than or equal to 1
    
    if n==0
        return "This method is not suited to use the matrix given"
    end

    m = standardform(m)
    
    x0 = zeros(size(A)[2]) # guess
    k = 0 # number of iterates

    while k<maxk

        residue = m-A*x0
        n0 = norm(residue, n)

        if 0.5*n0^2 < tol # convergence test
            x0=SquareReconstruction(xk)
            return A, m, x0, k
        end

        global x1 = x0 + A'*residue

        if norm(m-A*x1, n)>n0
            return A, m, x0, NaN # (||m-Axk||) should be a nonincreasing sequence
        end

        x0 = x1
        k=k+1

    end

    return A, m, x1, Inf

end


function SteepestDescent(filename; tol=1.e-4, maxk=100)

    A, m = loadfile(filename)  # retrieving the matrices

    m = standardform(m)

    V = Array{Matrix{Float64}}(undef, maxk+1)
    
    xk = A'*rand(size(A)[1]) # guess
    V[1] = SquareReconstruction(xk)
    k = 0 # number of iterates

    while k<maxk

        residue = A*xk-m

        n0 = 0.5*norm(residue)^2

        println(n0)

        if n0 < tol # convergence test
            xk=SquareReconstruction(xk)
            return A, m, xk, k, V
        end

        rk = A'*residue
        alphak=(norm(rk)/norm(A*rk))^2

        xk = xk - alphak*rk

        k=k+1

        V[k+1]=SquareReconstruction(xk)

    end

    xk=SquareReconstruction(xk)
    return A, m, xk, Inf, V

end