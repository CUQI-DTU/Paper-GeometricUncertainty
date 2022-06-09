function [res,f_vec] = lsqr_gentikh_weighted(A,D,b,u,alpha,lambda,options)
%This function solves the weighted generalized tikhonov regularized 
%least squares problem

%           min_x 1/2 ||(Ax-b)||_((lambda I + u^T u)^(-1))+alpha/2 ||Dx||_2^2
%           s.t.  x>=0 (optional)

%using the projected accelerated gradient descent method (with non-negative projections).
%The inverse covariance matrix (lambda I + u^T u)^(-1), where u is the model error samples.
%is multiplied using the woodbury formula for efficiency (see C_inv_mult
%function).

%INPUT:
%A: Forward Operator (M x N)
%D: Gaussian Prior Matrix (N x N)
%b: Projection data (M x 1)
%u: Model error sample matrix (M x Ns)
%alpha: Regularization parameter
%lambda: noise precision
%options: struct containing the fields:
    
    %maxiters: Maxiumum number of iterations (def. 50)
    %x0: Initial guess (def. zero vector)
    %iid: Boolean variable denoting if prior structure is identity (def. 1)
    %nonneg: Boolean variable denoting if we want non-negativity (def. 1)
    %epsilon: stopping tolerance (def. 10^(-3))

%Output:
%res: minimizer
%f_vec: vector of function values for different iterations

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
    h = 10^(-3);
end

uTu = u'*u;

%Determine Lipschitz constant using power iteration
n_lip = 20;
x = randn(n,1);

for i = 1:n_lip
    x1 = A'*(C_inv_mult(A*x,u,uTu,lambda));
    x = x1/norm(x1);
end

if iid == 1
    L = x'*A'*(C_inv_mult(A*x,u,uTu,lambda)) + delta;
else
    L = x'*A'*(C_inv_mult(A*x,u,uTu,lambda)) + 8*delta/h^2;
end

%Initialize algorithm
xold = x0;
y = x0;
told = 1;
k = 1;
converged = 0;

%Compute initial objective value
if nargout>1
    f_vec = 1/2*(A*x0-b)'*C_inv_mult(A*x0-b,u,uTu,lambda)+alpha/2*norm(D*x0,2)^2;
end

%Precompute A^T (lambda^(-1) I + uu^T)^(-1) b for efficiency
ATb = A'*C_inv_mult(b,u,uTu,lambda);

while k<maxiters && ~converged
    k = k+1;
    
    %Compute Gradient
    grad = A'*C_inv_mult(A*x,u,uTu,lambda)-ATb + alpha*(D'*D)*y;
    
    %Take Gradient step
    x = y-1/L*grad;
    
    %Implement non-negativity
    if nonneg
        x = max(0,x(:));
    end
    
    %Update t and y
    tnew = (1+sqrt(1+4*told^2))/2;
    y = x + (told-1)/tnew*(x-xold);
    
    %Compute relative difference in of norms between the iterates to
    %determine convergence
    relDif = norm(x-xold,2)/norm(xold,2);
    if relDif<epsilon
        converged = 1;
        disp(['The algorithm stopped at iteration ' num2str(k) ' because the relative change was below the specified tolerance'])
    end
    
    told = tnew;
    xold = x;
    
    %Compute objective value
    if nargout>1
        fval = 1/2*(A*x-b)'*C_inv_mult(A*x-b,u,uTu,lambda)+alpha/2*norm(D*x,2)^2;
        f_vec = [f_vec fval];
    end
end

if k==maxiters
    disp('The algorithm stopped because the number of iterations reached the specified maximum')
end
res = x;
end