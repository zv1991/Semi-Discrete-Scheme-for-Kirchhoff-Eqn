function simp = SimpsRule (h, u)
  simp = (h / 3) * (u (1) + u (end) + 4 * sum (u (2:2:end - 1)) + ...
  2 * sum (u (3:2:end - 2)));
endfunction
