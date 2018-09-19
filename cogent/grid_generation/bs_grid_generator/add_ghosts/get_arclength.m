function arc_length = get_arclength(R,Z)
  assert(length(R)==length(Z));

  arc_length(1) = 0.;
  for i=2:length(R)
    arc_length(i) = arc_length(i-1) + sqrt( (R(i) - R(i-1))^2 + (Z(i) - Z(i-1))^2 );
  end

end
