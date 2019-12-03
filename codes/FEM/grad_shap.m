%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  GRAD_SHAP computes the derivatives of the linear shape functions
%
%     grad_N = grad_shap(x)
%
%     Input:  x ... points
%
%     Output: grad_N ... derivatives of shape functions (dim nx2)

function grad_N = grad_shap(x)

  n = size(x,1);

  % preallocate memory
  grad_N = zeros(n,2);

  % compute values of derivatives
  grad_N(:,1) = -1/2*ones(n,1);
  grad_N(:,2) = 1/2*ones(n,1);

return