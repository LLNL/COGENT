function smooth_eqdsk(input_file, s)

  % If a negative value is passed for the smoothing parameter,
  % the algorithm will try to select a value automatically

  % Open the files
  fd_in  = fopen(input_file,'r');

  output_file = strcat(input_file,'_smoothed');
  fd_out = fopen(output_file,'w');

  dct_file = strcat(input_file,'_dct');
  fd_dct = fopen(dct_file,'w');

  % Open files
  fd_in  = fopen(input_file,'r');
  fd_out = fopen(output_file,'w');
  fd_dct = fopen(dct_file,'w');

  % Read and write the first line.  All we need is the problem size
  tline = fgets(fd_in);
  [header,count,errmsg,nextind] = sscanf(tline, '%s', 6);
  n = sscanf(tline(nextind:end), '%d', 2);
  nw = n(1);
  nh = n(2);
  fwrite(fd_out, tline);

  % Read and write the second and third lines containing the geometry data
  % we want
  tline = fgets(fd_in);
  geom1 = sscanf(tline, '%g', 5);
  fwrite(fd_out, tline);
  tline = fgets(fd_in);
  geom2 = sscanf(tline, '%g', 5);
  fwrite(fd_out, tline);

  % Read and write some lines we don't care about
  for i=1:2 + 4*ceil(nw/5)
    tline = fgets(fd_in);
    fwrite(fd_out, tline);
  end

  % Read the flux
  psi_orig = fscanf(fd_in, '%g', [nw nh])';

  % Read more data including safety factor,
  % LCFS and Limiter boundary
  % Seems that newer versions of g-file may have 
  % different formal, so use with caution
  read_extra_data = true
  if read_extra_data == true 
    % Read the safety factor
    safety_factor = fscanf(fd_in, '%g', nw);

    % Read the boundary (LCFS) and limiter points
    n = fscanf(fd_in, '%d', 2);
    nbdry = n(1);
    nlim = n(2);
    if nbdry>0
        boundary = fscanf(fd_in, '%g', [2 nbdry])';
        xb=boundary(:,1);
        yb=boundary(:,2);
    end
    if nlim>0
       limiter = fscanf(fd_in, '%g', [2 nlim])';
       xlim=limiter(:,1);
       ylim=limiter(:,2);
    end
  end
  
  % Get the geometry data
  rdim = geom1(1);
  zdim = geom1(2);
  rcentr = geom1(3);
  rleft = geom1(4);
  zmid = geom1(5);
  rmaxis = geom2(1);
  zmaxis = geom2(2);
  rmin = rleft;
  rmax = rmin + rdim;
  zmin = zmid - 0.5*zdim;
  zmax = zmid + 0.5*zdim;

  % Create linear space
  r = linspace(rmin, rmax, nw);
  z = linspace(zmin, zmax, nh);

  % Flip single-null geometry. 
  % In order to turn the upper SN to lower SN
  % configuration, we reflect psi w/r/t to z. 
  % We also need to change its sign  
  % to maintan proper direction of Bp 
  % e.g., toward the Xpoint.  
  use_reflected = false;
  if use_reflected 
    psi_lower_SN=zeros(nh,nw);
    for jj = 1:nh
        for ii = 1:nw
          psi_lower_SN(jj,ii) = -psi_orig(nh+1-jj,ii);
        end
    end
  else
    psi_lower_SN = psi_orig;  
  end  

  % set temporary psi to lower SN flux function  
  psi = psi_lower_SN;  

  % Add external point sources (e.g., to suppress upper Xpoint)
  use_external_sources = false;
  if use_external_sources 
    psi = psi_with_sources(psi_lower_SN,[rmin, rmax, zmin, zmax],[nw,nh]);
  end

  % Create flux function least-square fit 
  % e.g., to have a smoother psi in the divertor leg region
  use_fit = false;
  if use_fit 
    psi_fit = fitted_psi(psi,[rmin, rmax, zmin, zmax],[nw,nh]);
    psi = psi_fit;
  end  

  
  % Apply Garcia's smoothing algorithm 
  if s >= 0.
    [psi_smoothed psi_dct] = smoothn_mod(psi, s);
  else
    [psi_smoothed psi_dct] = smoothn_mod(psi);
  end

  % Write the smoothed flux
  fprintf(fd_out, '%16.9e%16.9e%16.9e%16.9e%16.9e\n', psi_smoothed');

  % Find the X point, given an initial guess (that might need to be adjusted)
  Xpt = [0.75*rmin+0.25*rmax, 0.75*zmin+0.25*zmax];
  Xpt = fminsearch(@(x) poloidal_Bmag(x, psi_dct', [rmin, rmax, zmin, zmax]), Xpt);
  fprintf('The X point is at (R,Z) = (%f,%f)\n', Xpt(1), Xpt(2));

  % Find the O point, given an initial guess (that might need to be adjusted)
  Opt = [rmaxis, zmaxis];
  Opt = fminsearch(@(x) poloidal_Bmag(x, psi_dct', [rmin, rmax, zmin, zmax]), Opt);
  fprintf('The O point is at (R,Z) = (%f,%f)\n', Opt(1), Opt(2));

  % Create normalized poloidal flux function
  psi_smoothed_Xpt = dct_derivs(Xpt(1), Xpt(2), psi_dct', [rmin, rmax, zmin, zmax],0);
  psi_smoothed_Opt = dct_derivs(Opt(1), Opt(2), psi_dct', [rmin, rmax, zmin, zmax],0);
  psi_smoothed_norm = (psi_smoothed - psi_smoothed_Opt)/(psi_smoothed_Xpt - psi_smoothed_Opt);

  % Write the dct coefficients
  fprintf(fd_dct, '%16.8f%16.8f%16.8f%16.8f%8d%8d\n', rmin, rmax, zmin, zmax, nw, nh);
  fprintf(fd_dct, '%16.8f%16.8f%16.8f%16.8f\n', rmaxis, zmaxis, Xpt(1), Xpt(2));
  fprintf(fd_dct, '%20.13e %20.13e %20.13e %20.13e %20.13e\n', psi_dct');


  % Make plots

  % Set selected values of the normalized flux function
  v = [1.15 1.1 1.05 1.0 0.9];

  figure(1);
  h1 = contour(r,z,psi_orig,80);
  if read_extra_data == true
    hold on
    plot(xlim,ylim, '-k', 'LineWidth',2)
    hold off
  end 
  cmin = min(psi_orig(:));
  cmax = max(psi_orig(:));
  caxis([cmin cmax]);
  colorbar
  title('Original flux function')
  axis image;

  figure(2);
  h2 = contour(r,z,psi_smoothed_norm,80);
  hold on
  contour(r,z,psi_smoothed_norm,v,'-k','LineWidth',1);
  hold off
  cmin = min(psi_smoothed_norm(:));
  cmax = max(psi_smoothed_norm(:));
  caxis([cmin cmax]);
  colorbar
  title('Normalized smoothed flux function')
  axis image;
 
  figure(3);
  func = (psi_smoothed - psi_lower_SN)/psi_smoothed_Xpt;
  h3 = pcolor(r,z,func);
  h3.FaceColor = 'interp';
  set(h3, 'EdgeColor', 'none');
  hold on
  contour(r,z,psi_smoothed_norm,v,'-k','LineWidth',1,'ShowText','on')
  hold off
  cmin=-0.1;
  cmax=0.1;
  caxis([cmin cmax]);
  colormap('turbo');
  colorbar
  title('Smoothing error estimate (psiSmoothed - psiLowerSN)/psiSmoothedXpt')
  axis image;


  % Read and write the rest of the file
  tline = fgets(fd_in);
  while ischar(tline)
    fwrite(fd_out, tline);
    tline = fgets(fd_in);
  end      

  % Close the files
  fclose(fd_in);
  fclose(fd_out);
  fclose(fd_dct);

end

function psi_fit = fitted_psi(psi, limits, Npts)

  % Create RBF interpolation by using point sources;
  % the idea here is that we can truncate the point source region
  % (e.g., zSmin) above the divertor legs, and then rely on smooth
  % (extrapolated) psi in the divertor leg region. We can do this
  % in case we need longer divertor legs that are supported 
  % by the original flux function.
  % NB: meaningful comparission between fitted 
  % and the GEQDSK data as well as modified GEQDSK file
  % is only possibe if extrapolated space
  % is identical to the original space 
  rmin = limits(1);
  rmax = limits(2); 
  zmin = limits(3);
  zmax = limits(4);
  nw = Npts(1);
  nh = Npts(2);
  r = linspace(rmin, rmax, nw);
  z = linspace(zmin, zmax, nh);

  % Set locations of point sources
  rSmin = limits(1);
  rSmax = limits(2);
  zSmin = limits(3);
  zSmax = limits(4);
  r0 = 0.5; %RBF is In/(r+r0)
  nSw = 20; %number of RBF's in r-direction
  nSh = 20; %number of RBF's in z-direction
  nS = nSw * nSh;
  zcoords = linspace(zSmin, zSmax, nSh);
  rcoords = linspace(rSmin, rSmax, nSw);
  nindex  = linspace(1,nS,nS);
  for i = 1:nSw
    for j = 1:nSh
     rS(j+(i-1)*nSh) = rcoords(i);
     zS(j+(i-1)*nSh) = zcoords(j);
    end
  end

  % Define point-source RBF's and weights of the original data 
  for nn = 1:nS
   for jj = 1:nh
    for ii = 1:nw
      psiPnt(jj,ii,nn) = 1.0/(sqrt((r(ii)-rS(nn))^2 + (z(jj)-zS(nn))^2) + r0); 
      weight(jj,ii) = 1.0;
      if ((z(jj) < zSmin) | (z(jj)>zSmax) | (r(ii)<rSmin) | (r(ii)>rSmax))
        weight(jj,ii) =0.0;
      end
    end
   end
  end 
  

  % Create least-square-fit matrix
  for n = 1:nS
   for k = 1:nS
     C(n,k) = sum(sum(weight.*psiPnt(:,:,n).*psiPnt(:,:,k)))/nw/nh;
   end
  end

  %Create the RHS of the linear system
  for n = 1:nS
    B(n) = sum(sum(weight.*psi.*psiPnt(:,:,n)))/nw/nh;
  end

  %Get point-source currents
  %opts.SYM = true;
  I = linsolve(C,B');

  %Define flux function fix on the extrapolated space
  psi_fit=zeros(nh,nw);
  for nn = 1:nS
   for jj = 1:nh
    for ii = 1:nw
      psi_fit(jj,ii) = psi_fit(jj,ii) + I(nn)/(sqrt((r(ii)-rS(nn))^2 + (z(jj)-zS(nn))^2) + r0);
    end
   end
  end

end

function psiWithSources = psi_with_sources(psi_lower_SN, limits, Npts)

   %Add extra point sources if needed
   %Format:
   % nSextr = n;
   % Iextr(1)=xxx; zSextr(1)=xxx; rSextr(1)=xxx; r0extr(1)=xxx;
   % Iextr(2)=xxx; zSextr(2)=xxx; rSextr(2)=xxx; r0extr(2)=xxx;
   % ...
   % Iextr(n)=xxx; zSextr(n)=xxx; rSextr(n)=xxx; r0extr(n)=xxx;

   nSextr=1;
   Iextr(1)=-0.005; zSextr(1)=1.5; rSextr(1)=1.2; r0extr(1)=0.1;

   psiWithSources = psi_lower_SN;

   rmin = limits(1);
   rmax = limits(2); 
   zmin = limits(3);
   zmax = limits(4);
   nw = Npts(1);
   nh = Npts(2);
   r = linspace(rmin, rmax, nw);
   z = linspace(zmin, zmax, nh);

   for nn = 1:nSextr
    for jj = 1:nh
     for ii = 1:nw
      psiWithSources(jj,ii) = psiWithSources(jj,ii) + Iextr(nn)/(sqrt((r(ii)-rSextr(nn))^2 + (z(jj)-zSextr(nn))^2) + r0extr(nn));
     end
    end
   end

end

function Bmag = poloidal_Bmag(x, coefs, limits)

  BR = -dct_derivs(x(1), x(2), coefs, limits, 2) / x(1);
  BZ =  dct_derivs(x(1), x(2), coefs, limits, 1) / x(1);

  Bmag = sqrt(BR^2 + BZ^2);

end

function sum = dct_derivs(u, v, f, limits, deriv)

  [N, M] = size(f);

  lambda = ones(max(N,M));
  lambda(1) = 1./sqrt(2.);

  Rmin = limits(1);
  Rmax = limits(2);
  Zmin = limits(3);
  Zmax = limits(4);

  scaled_u = (u-Rmin)*(N-1)/(Rmax-Rmin);
  scaled_v = (v-Zmin)*(M-1)/(Zmax-Zmin);

  sum = 0.;

  for j=1:M
    for i=1:N

      ufac1 = (i-1) * pi / N;
      vfac1 = (j-1) * pi / M;

      ufac2 = ufac1 * (scaled_u + 0.5);
      vfac2 = vfac1 * (scaled_v + 0.5);

      if deriv == 0
        sum = sum + lambda(i) * lambda(j) * f(i,j) * cos(ufac2) * cos(vfac2);
      elseif deriv == 1
        sum = sum - lambda(i) * lambda(j) * f(i,j) * sin(ufac2) * cos(vfac2) * ufac1;
      elseif deriv == 2
        sum = sum - lambda(i) * lambda(j) * f(i,j) * cos(ufac2) * sin(vfac2) * vfac1;
      elseif deriv == 3
        sum = sum + lambda(i) * lambda(j) * f(i,j) * sin(ufac2) * sin(vfac2) * ufac1 * vfac1;
      end
    end
  end

  sum = sum * 2. / sqrt(N*M);

  if deriv==1
    sum = sum * (N-1) / (Rmax-Rmin);
  elseif deriv==2
    sum = sum * (M-1) / (Zmax-Zmin);
  end

end

%{
  // Code to compute BZ, if needed for debugging  
  BZ=zeros(nh,nw);
  for jj = 2:nh-1
    for ii = 2:nw-1
      BZ(jj,ii) = (psi(jj,ii+1)-psi(jj,ii-1))/(r(ii+1)-r(ii-1))/r(ii);
    end
  end
%}