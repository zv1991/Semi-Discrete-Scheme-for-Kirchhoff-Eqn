## Author: Zurab Vashakidze <zurab.vashakidze@gmail.com>
## Created: 2023-07-03

## Legendre Polynomials
function [P_m] = LegendrePoly (n, x)
  for i = 0:n
    switch i
      case 0
        P_0 = 1;
        P_m = P_0;
      case 1
        P_1 = x;
        P_m = P_1;
      otherwise
        P_m = ((2 * i - 1) * x .* P_1 - (i - 1) * P_0) / i;
        P_0 = P_1;
        P_1 = P_m;
    endswitch
  endfor
endfunction
