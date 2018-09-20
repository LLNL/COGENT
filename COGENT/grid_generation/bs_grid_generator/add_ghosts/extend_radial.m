function [R_out,Z_out] = extend_radial(R_in,Z_in,n_extend,degree)

  n = length(R_in);

  % Make an arc length coordinate

  t_radial = get_arclength(R_in,Z_in);        

  npts = degree+1;

  R = R_in(n-npts+1:n);
  Z = Z_in(n-npts+1:n);
  t = t_radial(n-npts+1:n)';

  pp_R = polyfit(t, R, degree);
  pp_Z = polyfit(t, Z, degree);

  % Extrapolate the arc length coordinate

  t_uniform = linspace(0.,1.,npts)';        
  pp_arc = polyfit(t_uniform, t, degree);
  t_extend = 1. + linspace(1,n_extend,n_extend)' / (npts-1);
  t_arc_extend = polyval(pp_arc,t_extend);

  R_out  = polyval(pp_R,t_arc_extend);
  Z_out  = polyval(pp_Z,t_arc_extend);

end
