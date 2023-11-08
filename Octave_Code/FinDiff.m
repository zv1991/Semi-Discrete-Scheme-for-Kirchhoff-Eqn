function diff1_u = FinDiff (h, u)
  cfd = [1, -8, 0, 8, -1];
  ## the central differences %
  ffd = [-25, 48, -36, 16, -3];
  ## the forward differences %
  bfd = -ffd;
  ## the backward differences %
  len = length (u);
  ## determines length %
  diff1_u = zeros (1, len);
  ## preallocates vector %
  for j = 1:2
    diff1_u (j) = dot (ffd, u (j:j + 4)) / (12 * h);
    diff1_u (end + 1 - j) = dot (bfd, u (end + 1 - j:-1:end - 3 - j)) / ...
    (12 * h);
  endfor
  clear ('j')
  for j = 3:len - 2
    diff1_u (j) = dot (cfd, u (j - 2:j + 2)) / (12 * h);
  endfor
  clear ('j')
endfunction
