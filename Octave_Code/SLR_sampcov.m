## S_xx and S_xy as the sample variance and sample covariance, respectively

function S = SLR_sampcov (x, y)
  S = sum((x - mean(x) * ones(length(x),1)) .* (y - mean(y) * ...
  ones(length(y),1)));
endfunction
