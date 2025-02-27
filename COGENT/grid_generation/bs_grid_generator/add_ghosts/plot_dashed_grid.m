function plot_dashed_grid(data_R, data_Z, rgb, lw, radial, poloidal)

hold on;

[n_radial,n_poloidal] = size(data_R);

if radial

   t_poloidal_orig  = linspace(0.,1.,n_poloidal);
   tt_poloidal_fine = linspace(0.,1.,10*(n_poloidal));

   for i = 1:n_radial
      R_interp = spline(t_poloidal_orig, data_R(i,:), tt_poloidal_fine);
      Z_interp = spline(t_poloidal_orig, data_Z(i,:), tt_poloidal_fine);
      plot(R_interp,Z_interp,'Color',rgb/255,'LineStyle','--','LineWidth',lw);
   end
end

if poloidal

   t_radial_orig  = linspace(0.,1.,n_radial);
   tt_radial_fine = linspace(0.,1.,10*(n_radial));

   for j = 1:n_poloidal
      R_interp = spline(t_radial_orig, data_R(:,j), tt_radial_fine);
      Z_interp = spline(t_radial_orig, data_Z(:,j), tt_radial_fine);
      plot(R_interp,Z_interp,'Color',rgb/255,'LineStyle','--','LineWidth',lw);
   end
end

hold off;

end
