function dist = average_distance(R,Z)

  N = numel(R);

  dist = 0.;
  for i = 1:N-1
     dist = dist + ((R(i) - R(i+1))^2 + (Z(i) - Z(i+1))^2);
  end
  dist = sqrt(dist) / N;

end