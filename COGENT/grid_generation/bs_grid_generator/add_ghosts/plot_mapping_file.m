function plot_mapping_file(file_name)

   fd=fopen(file_name,'r');

   % Read LCORE
   bname = fscanf(fd,'%s',1);
   [nr_lcore, nr_lcore_extend, np_lcore, np_lcore_extend]=deal_array(fscanf(fd,'%d',4));
   [R_lcore,Z_lcore] = read_block(fd, nr_lcore, nr_lcore_extend, np_lcore, np_lcore_extend);

   fprintf('Block %s: nradial = %d, nradial_extend = %d, npoloidal = %d, npoloidal_extend = %d\n', bname, nr_lcore, nr_lcore_extend, np_lcore, np_lcore_extend);

   % Read RCORE
   bname = fscanf(fd,'%s',1);
   [nr_rcore, nr_rcore_extend, np_rcore, np_rcore_extend]=deal_array(fscanf(fd,'%d',4));
   [R_rcore,Z_rcore] = read_block(fd, nr_rcore, nr_rcore_extend, np_rcore, np_rcore_extend);

   fprintf('Block %s: nradial = %d, nradial_extend = %d, npoloidal = %d, npoloidal_extend = %d\n', bname, nr_rcore, nr_rcore_extend, np_rcore, np_rcore_extend);

   % Read LCSOL
   bname = fscanf(fd,'%s',1);
   [nr_lcsol, nr_lcsol_extend, np_lcsol, np_lcsol_extend]=deal_array(fscanf(fd,'%d',4));
   [R_lcsol,Z_lcsol] = read_block(fd, nr_lcsol, nr_lcsol_extend, np_lcsol, np_lcsol_extend);

   fprintf('Block %s: nradial = %d, nradial_extend = %d, npoloidal = %d, npoloidal_extend = %d\n', bname, nr_lcsol, nr_lcsol_extend, np_lcsol, np_lcsol_extend);

   % Read RCSOL
   bname = fscanf(fd,'%s',1);
   [nr_rcsol, nr_rcsol_extend, np_rcsol, np_rcsol_extend]=deal_array(fscanf(fd,'%d',4));
   [R_rcsol,Z_rcsol] = read_block(fd, nr_rcsol, nr_rcsol_extend, np_rcsol, np_rcsol_extend);

   fprintf('Block %s: nradial = %d, nradial_extend = %d, npoloidal = %d, npoloidal_extend = %d\n', bname, nr_rcsol, nr_rcsol_extend, np_rcsol, np_rcsol_extend);

   % Read LSOL
   bname = fscanf(fd,'%s',1);
   [nr_lsol, nr_lsol_extend, np_lsol, np_lsol_extend]=deal_array(fscanf(fd,'%d',4));
   [R_lsol,Z_lsol] = read_block(fd, nr_lsol, nr_lsol_extend, np_lsol, np_lsol_extend);

   fprintf('Block %s: nradial = %d, nradial_extend = %d, npoloidal = %d, npoloidal_extend = %d\n', bname, nr_lsol, nr_lsol_extend, np_lsol, np_lsol_extend);

   % Read RSOL
   bname = fscanf(fd,'%s',1);
   [nr_rsol, nr_rsol_extend, np_rsol, np_rsol_extend]=deal_array(fscanf(fd,'%d',4));
   [R_rsol,Z_rsol] = read_block(fd, nr_rsol, nr_rsol_extend, np_rsol, np_rsol_extend);

   fprintf('Block %s: nradial = %d, nradial_extend = %d, npoloidal = %d, npoloidal_extend = %d\n', bname, nr_rsol, nr_rsol_extend, np_rsol, np_rsol_extend);

   % Read LPF
   bname = fscanf(fd,'%s',1);
   [nr_lpf, nr_lpf_extend, np_lpf, np_lpf_extend]=deal_array(fscanf(fd,'%d',4));
   [R_lpf,Z_lpf] = read_block(fd, nr_lpf, nr_lpf_extend, np_lpf, np_lpf_extend);

   fprintf('Block %s: nradial = %d, nradial_extend = %d, npoloidal = %d, npoloidal_extend = %d\n', bname, nr_lpf, nr_lpf_extend, np_lpf, np_lpf_extend);

   % Read RPF
   bname = fscanf(fd,'%s',1);
   [nr_rpf, nr_rpf_extend, np_rpf, np_rpf_extend]=deal_array(fscanf(fd,'%d',4));
   [R_rpf,Z_rpf] = read_block(fd, nr_rpf, nr_rpf_extend, np_rpf, np_rpf_extend);

   fprintf('Block %s: nradial = %d, nradial_extend = %d, npoloidal = %d, npoloidal_extend = %d\n', bname, nr_rpf, nr_rpf_extend, np_rpf, np_rpf_extend);

   % Read MCORE
   bname = fscanf(fd,'%s',1);
   [nr_mcore, nr_mcore_extend, np_mcore, np_mcore_extend]=deal_array(fscanf(fd,'%d',4));
   [R_mcore,Z_mcore] = read_block(fd, nr_mcore, nr_mcore_extend, np_mcore, np_mcore_extend);

   fprintf('Block %s: nradial = %d, nradial_extend = %d, npoloidal = %d, npoloidal_extend = %d\n', bname, nr_mcore, nr_mcore_extend, np_mcore, np_mcore_extend);

   % Read MCSOL
   bname = fscanf(fd,'%s',1);
   [nr_mcsol, nr_mcsol_extend, np_mcsol, np_mcsol_extend]=deal_array(fscanf(fd,'%d',4));
   [R_mcsol,Z_mcsol] = read_block(fd, nr_mcsol, nr_mcsol_extend, np_mcsol, np_mcsol_extend);

   fprintf('Block %s: nradial = %d, nradial_extend = %d, npoloidal = %d, npoloidal_extend = %d\n', bname, nr_mcsol, nr_mcsol_extend, np_mcsol, np_mcsol_extend);

   lw = .25;

   new_R_lcore = reshape(R_lcore,nr_lcore+2*nr_lcore_extend,np_lcore+2*np_lcore_extend);
   new_Z_lcore = reshape(Z_lcore,nr_lcore+2*nr_lcore_extend,np_lcore+2*np_lcore_extend);
   R_lcore_valid = new_R_lcore(nr_lcore_extend+1:nr_lcore_extend+nr_lcore,np_lcore_extend+1:np_lcore_extend+np_lcore);
   Z_lcore_valid = new_Z_lcore(nr_lcore_extend+1:nr_lcore_extend+nr_lcore,np_lcore_extend+1:np_lcore_extend+np_lcore);
   plot_grid(R_lcore_valid,Z_lcore_valid, [0 0 0], lw, true, true);

   new_R_mcore = reshape(R_mcore,nr_mcore+2*nr_mcore_extend,np_mcore+2*np_mcore_extend);
   new_Z_mcore = reshape(Z_mcore,nr_mcore+2*nr_mcore_extend,np_mcore+2*np_mcore_extend);
   R_mcore_valid = new_R_mcore(nr_mcore_extend+1:nr_mcore_extend+nr_mcore,np_mcore_extend+1:np_mcore_extend+np_mcore);
   Z_mcore_valid = new_Z_mcore(nr_mcore_extend+1:nr_mcore_extend+nr_mcore,np_mcore_extend+1:np_mcore_extend+np_mcore);
   plot_grid(R_mcore_valid,Z_mcore_valid, [0 0 0], lw, true, true);

   new_R_rcore = reshape(R_rcore,nr_rcore+2*nr_rcore_extend,np_rcore+2*np_rcore_extend);
   new_Z_rcore = reshape(Z_rcore,nr_rcore+2*nr_rcore_extend,np_rcore+2*np_rcore_extend);
   R_rcore_valid = new_R_rcore(nr_rcore_extend+1:nr_rcore_extend+nr_rcore,np_rcore_extend+1:np_rcore_extend+np_rcore);
   Z_rcore_valid = new_Z_rcore(nr_rcore_extend+1:nr_rcore_extend+nr_rcore,np_rcore_extend+1:np_rcore_extend+np_rcore);
   plot_grid(R_rcore_valid,Z_rcore_valid, [0 0 0], lw, true, true);

   new_R_lcsol = reshape(R_lcsol,nr_lcsol+2*nr_lcsol_extend,np_lcsol+2*np_lcsol_extend);
   new_Z_lcsol = reshape(Z_lcsol,nr_lcsol+2*nr_lcsol_extend,np_lcsol+2*np_lcsol_extend);
   R_lcsol_valid = new_R_lcsol(nr_lcsol_extend+1:nr_lcsol_extend+nr_lcsol,np_lcsol_extend+1:np_lcsol_extend+np_lcsol);
   Z_lcsol_valid = new_Z_lcsol(nr_lcsol_extend+1:nr_lcsol_extend+nr_lcsol,np_lcsol_extend+1:np_lcsol_extend+np_lcsol);
   plot_grid(R_lcsol_valid,Z_lcsol_valid, [0 0 0], lw, true, true);

   new_R_mcsol = reshape(R_mcsol,nr_mcsol+2*nr_mcsol_extend,np_mcsol+2*np_mcsol_extend);
   new_Z_mcsol = reshape(Z_mcsol,nr_mcsol+2*nr_mcsol_extend,np_mcsol+2*np_mcsol_extend);
   R_mcsol_valid = new_R_mcsol(nr_mcsol_extend+1:nr_mcsol_extend+nr_mcsol,np_mcsol_extend+1:np_mcsol_extend+np_mcsol);
   Z_mcsol_valid = new_Z_mcsol(nr_mcsol_extend+1:nr_mcsol_extend+nr_mcsol,np_mcsol_extend+1:np_mcsol_extend+np_mcsol);
   plot_grid(R_mcsol_valid,Z_mcsol_valid, [0 0 0], lw, true, true);

   new_R_rcsol = reshape(R_rcsol,nr_rcsol+2*nr_rcsol_extend,np_rcsol+2*np_rcsol_extend);
   new_Z_rcsol = reshape(Z_rcsol,nr_rcsol+2*nr_rcsol_extend,np_rcsol+2*np_rcsol_extend);
   R_rcsol_valid = new_R_rcsol(nr_rcsol_extend+1:nr_rcsol_extend+nr_rcsol,np_rcsol_extend+1:np_rcsol_extend+np_rcsol);
   Z_rcsol_valid = new_Z_rcsol(nr_rcsol_extend+1:nr_rcsol_extend+nr_rcsol,np_rcsol_extend+1:np_rcsol_extend+np_rcsol);
   plot_grid(R_rcsol_valid,Z_rcsol_valid, [0 0 0], lw, true, true);

   new_R_lsol = reshape(R_lsol,nr_lsol+2*nr_lsol_extend,np_lsol+2*np_lsol_extend);
   new_Z_lsol = reshape(Z_lsol,nr_lsol+2*nr_lsol_extend,np_lsol+2*np_lsol_extend);
   R_lsol_valid = new_R_lsol(nr_lsol_extend+1:nr_lsol_extend+nr_lsol,np_lsol_extend+1:np_lsol_extend+np_lsol);
   Z_lsol_valid = new_Z_lsol(nr_lsol_extend+1:nr_lsol_extend+nr_lsol,np_lsol_extend+1:np_lsol_extend+np_lsol);
   plot_grid(R_lsol_valid,Z_lsol_valid, [0 0 0], lw, true, true);

   new_R_rsol = reshape(R_rsol,nr_rsol+2*nr_rsol_extend,np_rsol+2*np_rsol_extend);
   new_Z_rsol = reshape(Z_rsol,nr_rsol+2*nr_rsol_extend,np_rsol+2*np_rsol_extend);
   R_rsol_valid = new_R_rsol(nr_rsol_extend+1:nr_rsol_extend+nr_rsol,np_rsol_extend+1:np_rsol_extend+np_rsol);
   Z_rsol_valid = new_Z_rsol(nr_rsol_extend+1:nr_rsol_extend+nr_rsol,np_rsol_extend+1:np_rsol_extend+np_rsol);
   plot_grid(R_rsol_valid,Z_rsol_valid, [0 0 0], lw, true, true);

   new_R_lpf = reshape(R_lpf,nr_lpf+2*nr_lpf_extend,np_lpf+2*np_lpf_extend);
   new_Z_lpf = reshape(Z_lpf,nr_lpf+2*nr_lpf_extend,np_lpf+2*np_lpf_extend);
   R_lpf_valid = new_R_lpf(nr_lpf_extend+1:nr_lpf_extend+nr_lpf,np_lpf_extend+1:np_lpf_extend+np_lpf);
   Z_lpf_valid = new_Z_lpf(nr_lpf_extend+1:nr_lpf_extend+nr_lpf,np_lpf_extend+1:np_lpf_extend+np_lpf);
   plot_grid(R_lpf_valid,Z_lpf_valid, [0 0 0], lw, true, true);

   new_R_rpf = reshape(R_rpf,nr_rpf+2*nr_rpf_extend,np_rpf+2*np_rpf_extend);
   new_Z_rpf = reshape(Z_rpf,nr_rpf+2*nr_rpf_extend,np_rpf+2*np_rpf_extend);
   R_rpf_valid = new_R_rpf(nr_rpf_extend+1:nr_rpf_extend+nr_rpf,np_rpf_extend+1:np_rpf_extend+np_rpf);
   Z_rpf_valid = new_Z_rpf(nr_rpf_extend+1:nr_rpf_extend+nr_rpf,np_rpf_extend+1:np_rpf_extend+np_rpf);
   plot_grid(R_rpf_valid,Z_rpf_valid, [0 0 0], lw, true, true);

   axis equal;
   xlim([.5 2.75]);
   ylim([-1.65 1.25]);
   xlabel('R (m)');
   ylabel('Z (m)');
end

function [R, Z] = read_block(fd, nradial, nradial_extend, npoloidal, npoloidal_extend)

   num_pts = (nradial + 2*nradial_extend) * (npoloidal + 2*npoloidal_extend);

   for i=1:num_pts
      R(i) = fscanf(fd,'%g', 1);
      Z(i) = fscanf(fd,'%g', 1);
      dummy = fscanf(fd,'%g', 1);
      dummy = fscanf(fd,'%g', 1);
   end

end

function varargout = deal_array(arr)
    s = numel(arr);
    n = nargout;

    if n > s
        error('Insufficient number of elements in array!');
    elseif n == 0
        return;
    end

    for i = 1:n
        varargout(i) = {arr(i)}; %#ok<AGROW>
    end
end

function plot_grid(data_R, data_Z, rgb, lw, radial, poloidal)

hold on;

[n_radial,n_poloidal] = size(data_R);

if radial

   t_poloidal_orig  = linspace(0.,1.,n_poloidal);
   tt_poloidal_fine = linspace(0.,1.,10*(n_poloidal));

   for i = 1:n_radial
      R_interp = spline(t_poloidal_orig, data_R(i,:), tt_poloidal_fine);
      Z_interp = spline(t_poloidal_orig, data_Z(i,:), tt_poloidal_fine);
      plot(R_interp,Z_interp,'Color',rgb/255,'LineWidth',lw);
   end
end

if poloidal

   t_radial_orig  = linspace(0.,1.,n_radial);
   tt_radial_fine = linspace(0.,1.,10*(n_radial));

   for j = 1:n_poloidal
      R_interp = spline(t_radial_orig, data_R(:,j), tt_radial_fine);
      Z_interp = spline(t_radial_orig, data_Z(:,j), tt_radial_fine);
      plot(R_interp,Z_interp,'Color',rgb/255,'LineWidth',lw);
   end
end

hold off;

end
