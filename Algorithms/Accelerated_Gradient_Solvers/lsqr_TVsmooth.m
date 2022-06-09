function [x,f_vec] = lsqr_TVsmooth(A,b,alpha,options)

%This function solves a smooth approximation to the TV regularized least squares problem

%   min_x 1/2||Ax-b||_2^2+alpha/2 TV_{rho}(x)
%   s.t. x>=0

%where the regularization term TV_{rho} is a smooth approximation to the TV
%of x. The problem is solved using accelerated projected gradient descent.

%INPUT:
%A: System matrix
%b: data
%alpha: regularization parameter
%options: struct with the fields

    %maxiters: Maximum number of iterations. (def. 50)
    %x0: initial guess. (def. vector of zeros)
    %nonneg: boolean parameter for nonnegativity. (def. 1)
    %rho: TV smoothness parameter. (def. 10^(-4))
    %h: Finite difference discretization parameter. (def. 10^(-3))
    %epsilon: Tolerance parameter. (def. 10^(-4))

%OUTPUT:
%x: reconstruction
%f_vec: vector of function values

[~,n] = size(A);
N = sqrt(n);

if nargin>3
    maxiters = options.maxiters;
    x0 = options.x0;
    nonneg = options.nonneg;
    rho = options.rho;
    h = options.h;
    epsilon = options.epsilon;
else
    maxiters = 50;
    x0 = zeros(n,1);
    nonneg = 1;
    rho = 10^(-6);
    h = 10^(-2);
    epsilon = 10^(-6);
end

%Compute Lipschitz constant
L = normest(A)^2 + 8*alpha/rho/h^2;

%Compute Finite difference matrix with Neumann b.c.
Dfd = spdiags([-ones(N-1,1),ones(N-1,1);0 1]/h,[0,1],N,N);
D = vertcat(kron(Dfd,speye(N)),kron(speye(N),Dfd));
% Compute "smooth" TV and its gradient
phi = @(y) sqrt(sum(reshape(y,[],2).^2,2)+rho^2);

%Initialize
k = 0;
x = x0;
y = x;
xold = x;
converged = 0;
told = 1;

Dx = D*x;
smooth_tv = sum(phi(Dx));

if nargout>1
    fval = 1/2*norm(A*x-b,2)^2+alpha/2*smooth_tv;
    f_vec = fval;
end

while ~converged && k<maxiters
    k = k+1;
    
    %Compute TV-gradient
    Dy = D*y;
    grad_tv = D'*(Dy.*repmat(1./phi(Dy),2,1));
    
    %Take gradient step
    res = A*y-b;
    grad = A'*res+alpha*grad_tv;
    
    x = y - (1/L)*grad;
    
    %Do a projection step
    if nonneg
        x = max(0,x);
    end 
    
    %Do momentum step
    t = (1 + sqrt(1+4*told^2))/2; 
    y = x + (told-1)/t*(x-xold);
    
    if norm(x-xold)/norm(xold)<epsilon
        disp(['The algorithm stopped at iteration ' num2str(k) ' because the relative change in objective is below the set tolerance'])
        converged = 1;
    end
    
    xold = x;
    told = t;
    
    if nargout>1
        smooth_tv = sum(phi(Dx));
        fval = 1/2*norm(A*x-b,2)^2+alpha/2*smooth_tv;
        f_vec = [f_vec,fval];
    end
end
if k==maxiters
    disp('The algorithm stopped because the number of iterations reached the set maximum')
end
end