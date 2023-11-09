## Simple Linear Regression

function [beta1, beta0, y2, res, min_res, max_res, r_squared] = SLR (x, y)
  S_xy      = SLR_sampcov (x, y);
  S_xx      = SLR_sampcov (x, x);
  S_yy      = SLR_sampcov (y, y);
  beta1     = S_xy / S_xx;
  beta0     = mean(y) - beta1*mean(x);
  y2        = beta0 + x*beta1;
  res       = y - y2;
  min_res   = min(abs(res));
  max_res   = max(abs(res));
  r_squared = (S_xy * S_xy) / (S_xx * S_yy);
endfunction
