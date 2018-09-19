function d = deriv(f)

  for i = 2:length(f)
    d(i) = f(i) - f(i-1);
  end

  d(1) = d(2);

end
