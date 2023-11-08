function sol = tridiagonal_algotirhm (b, phi)
  ## determines len
  len = length (b);
  ## preallocates solution vector
  sol = zeros (1, len);
  ## forward elimination
  for i = 2:len
    temp = 1 / b (i - 1);
    b (i) = b (i) - temp;
    phi (i) = phi (i) - phi (i - 1) * temp;
  endfor
  clear ('i')
  ## backward substitution
  sol (len) = phi (len) / b (len);
  for i = (len - 1):(-1):1
    sol (i) = (phi (i) - sol (i + 1)) / b (i);
  endfor
  clear ('i')
endfunction

