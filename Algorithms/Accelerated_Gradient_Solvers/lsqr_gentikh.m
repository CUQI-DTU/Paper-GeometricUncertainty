function [res,f_vec] = lsqr_gentikh(A,D,b,u,lambda,delta,options)
%This function solves the Generalized Tikhonov Minimization Problem

%           min_x lambda/2 ||Ax-b||_2^2+delta/2 ||Dx-u||_2^2
%           s.t.  x\in C

%using the FISTA method. (C correponds to either no constraints or
%non-negativity (controlled in options))

%INPUT:
%A: Forward Operator (M x N "matrix")
%D: Gaussian Prior Precision matrix (N x N "matrix)
%b: Vectorized sinogram (M x 1 vector)
%u: Prior mean shift    (N x 1 vector)
%lambda: noise precision parameter
%delta: Gaussian prior parameter
%options: struct containing the fields:
    
    %maxiters: Maxiumum number of iterations (default 50)
    %x0: Initial guess (default zero vector)
    %iid: Boolean variable denoting if prior structure is identity (default 1)
    %nonneg: Boolean variable denoting if we want non-negativity (default 1)
    %epsilon: Heuristic stopping tolerance (def. 10^(-8))
    %h: precision matrix discretization length (default 1)

%OUTPUT:
%res: minimizer
%f_vec: vector of objective values for different iterations

n = size(D,1);

if nargin>6
    x0 = options.x0;
    maxiters = options.maxiters;
    nonneg = options.nonneg;
    epsilon = options.epsilon;
    iid = options.iid;
    h = options.h;
else
    x0 = zeros(n,1);
    maxiters = 50;
    nonneg = 1;
    epsilon = 10^(-8);
    iid = 1;
    h = 1;
end

%Determine Lipschitz constant using power iteration

if iid == 1
    L = lambda*normest(A)^2 + delta;
else
    L = lambda*normest(A)^2 + 8*delta/h^2;
end

%Initialize algorithm
xold = x0;
y = x0;
told = 1;
k = 1;
converged = 0;

%Compute initial objective value
if nargout>1
    f_vec = lambda/2*norm(A*x0-b,2)^2+delta/2*norm(x0-u,2)^2;
end


while k<maxiters && ~converged
    k = k+1;
    
    %Evaluate gradient
    grad = lambda*A'*(A*y-b)+delta*D'*(D*y-u);
    
    %Take gradient step
    x = y-1/L*grad;
    
    %Nonnegativity projection
    if nonneg == 1
        x = max(0,x(:));
    end
    
    %Update t and y including momentum term
    tnew = (1+sqrt(1+4*told^2))/2;
    y = x + (told-1)/tnew*(x-xold);
    
    told = tnew;
    xold = x;
    
    %Compute objective value
    if nargout>1
        fval = lambda/2*norm(A*x-b,2)^2+delta/2*norm(D*x-u,2)^2;
        f_vec = [f_vec fval];
    end
end

res = x;
end