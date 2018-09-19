function [interp_R, interp_Z, nr, nr_extend, np, np_extend] = add_block_ghosts(mapping_file, block)

  %%%% Set options %%%%%%%%%%%%%%%%%

  plot_valid_grid = true;
  plot_extrapolated_grid = false;
  plot_interpolated_grid = true;

  nr_extend = 4;
  np_extend = 4;

  rgb;
  if strcmp(block,'lcore')
    color_map = Red;
  elseif strcmp(block,'rcore')
    color_map = Green;      
  elseif strcmp(block,'mcore')
    color_map = Cyan;
  elseif strcmp(block,'lcsol')
    color_map = Tan;
  elseif strcmp(block,'rcsol')
    color_map = Blue;      
  elseif strcmp(block,'mcsol')
    color_map = Violet_Red;
  elseif strcmp(block,'lsol')
    color_map = Light_Blue;
  elseif strcmp(block,'rsol')
    color_map = Goldenrod;      
  elseif strcmp(block,'lpf')
    color_map = Orange;
  elseif strcmp(block,'rpf')
    color_map = Purple;  % purple      
  end

  r_extrap_degree = 1;
  p_extrap_degree = 2;

  %%%% End options %%%%%%%%%%%%%%%%%%

  fd=fopen(mapping_file,'r');

  % Find the block
  found_block = false;
  for block_num=1:10

    bname = fscanf(fd,'%s',1);
    [b_nr, b_nr_extend, b_np, b_np_extend]=deal_array(fscanf(fd,'%d',4));
    [data_R,data_Z] = read_valid_block(fd, b_nr, b_nr_extend, b_np, b_np_extend);

    if strcmp(block,bname)
      found_block = true;
      break;
    else
      clear data_R data_Z b_nr b_nr_extend b_np b_np_extend;
    end
  end

  fclose(fd);

  if ~found_block
    fprintf('Block %s not found in file %s\n', block, mapping_file);
    return;
  end

  fprintf('Adding ghost points to %s block\n', block);

  [nr np] = size(data_R);

  figure(1);
  axis ([0.7 2.5 -1.7 1.5]);
  xlabel('R');
  ylabel('Z','Rotation',0.);
  axis equal;

  if plot_valid_grid
    plot_grid(data_R, data_Z, color_map, 2., true, true);
  end

  % Extrapolate the valid grid poloidally

  pol_extrap_R = zeros(nr, np+2*np_extend);
  pol_extrap_Z = zeros(nr, np+2*np_extend);

  for i=1:nr

    % Copy interior
    for j=1:np
      pol_extrap_R(i,j+np_extend)  = data_R(i,j);
      pol_extrap_Z(i,j+np_extend)  = data_Z(i,j);
    end 

    if np_extend > 0

      % Extrapolate at low boundary
      [R_lo_extension,Z_lo_extension] = ... 
          extend_poloidal(fliplr(data_R(i,:)),fliplr(data_Z(i,:)),np_extend,p_extrap_degree);
      for j=1:np_extend
        pol_extrap_R(i,np_extend-j+1) = R_lo_extension(j);
        pol_extrap_Z(i,np_extend-j+1) = Z_lo_extension(j);
      end

      % Extrapolate at high boundary
      [R_hi_extension,Z_hi_extension] = extend_poloidal(data_R(i,:),data_Z(i,:),np_extend,p_extrap_degree);
      for j=1:np_extend
        pol_extrap_R(i,j+np+np_extend) = R_hi_extension(j);
        pol_extrap_Z(i,j+np+np_extend) = Z_hi_extension(j);
      end

    end
  end

  % Extrapolate the poloidally extrapolted grid radially

  extrap_R = zeros(nr+2*nr_extend, np+2*np_extend);
  extrap_Z = zeros(nr+2*nr_extend, np+2*np_extend);

  for j=1:np + 2*np_extend

    % Copy interior
    for i=1:nr
        extrap_R(i+nr_extend,j) = pol_extrap_R(i,j);
        extrap_Z(i+nr_extend,j) = pol_extrap_Z(i,j);
    end 

    if nr_extend > 0

      % Extrapolate at low boundary
      [R_lo_extension,Z_lo_extension] = ...
            extend_radial(flipud(pol_extrap_R(:,j)),flipud(pol_extrap_Z(:,j)),nr_extend,r_extrap_degree);
      for i=1:nr_extend
        extrap_R(nr_extend-i+1,j) = R_lo_extension(i);
        extrap_Z(nr_extend-i+1,j) = Z_lo_extension(i);
      end

      % Extrapolate at high boundary
      [R_hi_extension,Z_hi_extension] = extend_radial(pol_extrap_R(:,j),pol_extrap_Z(:,j),nr_extend,r_extrap_degree);
      for i=1:nr_extend
        extrap_R(nr + nr_extend + i,j)  = R_hi_extension(i);
        extrap_Z(nr + nr_extend + i,j)  = Z_hi_extension(i);
      end
    end      

  end

  [m0 m1] = size(extrap_R);

  npts=0;
  for j=1+np_extend:m1-np_extend
    for i=1+nr_extend:m0-nr_extend
      npts = npts+1;
      xi0(npts) = i;
      xi1(npts) = j;
      node_R(npts) = data_R(i-nr_extend,j-np_extend);
      node_Z(npts) = data_Z(i-nr_extend,j-np_extend);
    end
  end

  for j=np_extend+1:m1-2*np_extend-1:m1-np_extend
    for i=1:nr_extend
      npts = npts+1;
      xi0(npts) = i;
      xi1(npts) = j;
      node_R(npts) = extrap_R(i,j);
      node_Z(npts) = extrap_Z(i,j);
    end
    for i=m0-nr_extend+1:m0
      npts = npts+1;
      xi0(npts) = i;
      xi1(npts) = j;
      node_R(npts) = extrap_R(i,j);
      node_Z(npts) = extrap_Z(i,j);
    end
  end

  for i=nr_extend+1:m0-2*nr_extend-1:m0-nr_extend
    for j=1:np_extend
      npts = npts+1;
      xi0(npts) = i;
      xi1(npts) = j;
      node_R(npts) = extrap_R(i,j);
      node_Z(npts) = extrap_Z(i,j);
    end
    for j=m1-np_extend+1:m1
      npts = npts+1;
      xi0(npts) = i;
      xi1(npts) = j;
      node_R(npts) = extrap_R(i,j);
      node_Z(npts) = extrap_Z(i,j);
    end
  end

  if strcmp(block,'mcore') || strcmp(block,'mcsol')

    if strcmp(block,'mcore')
      step = 8;
    end 
    if strcmp(block,'mcsol') 
      step = 1;
    end

    for j=np_extend+1+step:step:m1-np_extend-step
      for i=1:nr_extend
        npts = npts+1;
        xi0(npts) = i;
        xi1(npts) = j;
        node_R(npts) = extrap_R(i,j);
        node_Z(npts) = extrap_Z(i,j);
      end
      for i=m0-nr_extend+1:m0
        npts = npts+1;
        xi0(npts) = i;
        xi1(npts) = j;
        node_R(npts) = extrap_R(i,j);
        node_Z(npts) = extrap_Z(i,j);
      end
    end

  end

  average_dist = average_distance(xi0,xi1);
  % fprintf('   Average distance between interpolation points = %e\n', average_dist);
  eps = 0.5 * average_dist;
  % fprintf('   RBF parameter = %e\n', eps);

  rbf_obj_R = rbfcreate([xi0;xi1], node_R, 'RBFFunction', 'cubic', 'RBFConstant', eps);
  rbf_obj_Z = rbfcreate([xi0;xi1], node_Z, 'RBFFunction', 'cubic', 'RBFConstant', eps);

  for j=1:m1
    for i=1:m0
      interp_R(i,j) = rbfinterp([i;j],rbf_obj_R);
      interp_Z(i,j) = rbfinterp([i;j],rbf_obj_Z);
    end
  end    

  if plot_interpolated_grid
    %plot_grid(interp_R, interp_Z, Black, 1., true, true);
    plot_dashed_grid(interp_R, interp_Z, Black, 1., true, true);
  end

  if plot_extrapolated_grid
    %plot_dashed_grid(extrap_R, extrap_Z, color_map, 0.5, true, true);
    plot_grid(extrap_R, extrap_Z, color_map, 0.5, true, true);
  end

  return;


function [R, Z] = read_valid_block(fd, nr, nr_extend, np, np_extend)

  for j=1:np + 2*np_extend
    for i=1:nr + 2*nr_extend
      this_R = fscanf(fd,'%g', 1);
      this_Z = fscanf(fd,'%g', 1);
      dummy = fscanf(fd,'%g', 1);
      dummy = fscanf(fd,'%g', 1);

      if(i > nr_extend && i <= nr + nr_extend && j > np_extend && j <= np + np_extend)
        R(i-nr_extend,j-np_extend) = this_R;
        Z(i-nr_extend,j-np_extend) = this_Z;
       end      

    end
  end

  return;


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

  return;

