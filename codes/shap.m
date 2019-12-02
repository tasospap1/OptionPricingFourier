%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SHAP computes the values of the shape functions for linear elements
%
%     N = shap(x)
%
%     Input:  x ... points
%
%     Output: N ... shape functions (dim nx2)

function N = shap(x)

  n = size(x,1);

  % preallocate memory
  N = zeros(n,2);

  % compute function values
  N(:,1) = 1/2*(1-x);
  N(:,2) = 1/2*(1+x);

return
