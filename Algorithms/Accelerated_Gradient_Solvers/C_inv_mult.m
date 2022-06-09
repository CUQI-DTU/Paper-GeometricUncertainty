function res = C_inv_mult(x,u,uTu,lambda)
%This function computes the matrix vector multiplication
%           (lambda^(-1) I + uu^T)^(-1) x
%by utilizing the Woodbury formula:
%       lambda I-lambda^2 u(I+u^T U)^(-1) u^T
%which is neccessary for efficiently utilizing the model-discrepancy method

%INPUT:
%x: Input vector of length M
%u: Model discrepency sample matrix (M x N_s)
%uTu: Matrix-Matrix product u^T u (N_s x N_s)
%lambda: Noise precision parameter

%Output:
%res: Resulting matrix vector multiplication

[~,Ns] = size(u);

%We start the multiplication from the right side:
uT_x = u'*x;

%Solve (N_s x N_s) linear system directly
tmp = (eye(Ns)+lambda*uTu)\uT_x;

%Finish it
res = lambda*x - lambda^2*u*tmp;
end
