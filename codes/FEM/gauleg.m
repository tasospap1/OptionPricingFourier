%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  GAULEG computes Gauss quadrature points
%
%     [x w] = gauleg(n)
%
%     Input:  n ... number of gauss points
%
%     Output: x ... gauss points
%             w ... weights

function [x,w] = gauleg(n)

  % Initalize variables
  m = floor((n+1)/2);
  x  = zeros(n,1);
  w  = zeros(n,1);

  for i = 1:m

    % Initial guess of root (starting value)
    z = cos(pi*(i-1/4)/(n+1/2));

    delta = 1;
    while(delta > eps)

      p1 = 0;
      p2 = 1;

      for k = 0:(n-1)

        % Computing value of n-th Legendre polynomial at point z using the
        % recursion:
        %
        %   (j+1)*P_(j+1)(z) = (2*j+1)*z*P_(j)(z)-j*P_(j-1)(z)

        p3 = ((2*k+1)*z*p2-k*p1)/(k+1);

        % Computing value of first derivative of n-th Legendre polynomial
        % at point z using the recursion:
        %
        %   (1-z^2)*P'_(j)(z) = j*[z*P_(j)(z)-P_(j-1)(z)]

        dp = n*(z*p3-p2)/(z^2-1);
        p1 = p2;
        p2 = p3;

      end

      % Performing Newton update

      z_old = z;
      z = z_old-p3/dp;

      delta = abs(z-z_old);

    end

    % Compute Gauss points in [-1 1]
    x(i) = -z;
    x(n+1-i) = z;
    w(i) = 2/((1-z^2)*dp^2);
    w(n+1-i) = w(i);

  end

  return
