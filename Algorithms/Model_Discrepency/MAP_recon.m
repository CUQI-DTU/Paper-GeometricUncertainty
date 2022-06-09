function [x_MAP,fvec] = MAP_recon(A,b,alpha,x0,reg_term,maxiters,nonneg)

%This function computes the solution to the optimization problem
%   min_x 1/2 ||Ax-b||^2+alpha/2*R(x),
%where R(x) is a regularization term (tikh, gentikh or TV). There is also a
%possibility to add nonnegativity contraints.

%INPUT:
%A: SYSTEM matrix
%b: data
%alpha: regularization parameter
%x0: initial guess
%options: solver options (regularization term, maximum iterations,
%nonnegativity)

%OUTPUT:
%x_MAP: MAP solution

if strcmp(reg_term,'TV')
    %Solve smooth TV problem using FISTA
    options.x0 = x0;
    options.maxiters = maxiters;
    options.rho = 10^(-2);
    options.epsilon = 10^(-6);
    options.h = 10^(-3);
    options.nonneg = nonneg;
    
    [x_MAP,fvec] = lsqr_TVsmooth(A,b,alpha,options);
    
elseif strcmp(reg_term,'gentikh')
    N = sqrt(length(x0));
    %Build precision matrix  
    e_vec = ones(N,1);
    L_1D = spdiags([-e_vec 2*e_vec -e_vec], [-1 0 1],N,N);
    L = kron(speye(N),L_1D) + kron(L_1D,speye(N));
    D = chol(L);
    
    options.h = 1;
    options.maxiters = maxiters;
    options.iid = 0;
    options.x0 = x0;
    options.epsilon = 10^(-6);
    options.nonneg = nonneg;
    
    u = zeros(N^2,1);
    
    [x_MAP,fvec] = lsqr_gentikh(A,D,b,u,1,alpha,options);
elseif strcmp(reg_term,'tikh')
    N = sqrt(length(x0));
   %Build precision matrix
    D = speye(N^2);
    
    options.h = 1;
    options.maxiters = maxiters;
    options.iid = 0;
    options.x0 = x0;
    options.epsilon = 10^(-6);
    options.nonneg = nonneg;
    
    u = zeros(N^2,1);
    
    [x_MAP,fvec] = lsqr_gentikh(A,D,b,u,1,alpha,options);
end
end