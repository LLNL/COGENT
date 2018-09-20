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
  psi = fscanf(fd_in, '%g', [nw nh])';

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

  r = linspace(rmin, rmax, nw);
  z = linspace(zmin, zmax, nh);

  % Created extrapolated space
  % NB: meaningful comparission between fitted 
  % and the GEQDSK data as well as modified GEQDSK file
  % is only possibe if extrapolated space
  % is identical to the original space 
  rminExt = rmin;
  rmaxExt = rmax; 
  zminExt = zmin;
  zmaxExt = zmax;
  nwExt = nw;
  nhExt = nh;
  rExt = linspace(rminExt, rmaxExt, nwExt);
  zExt = linspace(zminExt, zmaxExt, nhExt);

  % Set locations of point sources
  zSmin = zmin;
  zSmax = zmax;
  rSmin = rmin;
  rSmax = rmax;
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
  psi_fit=zeros(nhExt,nwExt);
  for nn = 1:nS
   for jj = 1:nhExt
    for ii = 1:nwExt
      psi_fit(jj,ii) = psi_fit(jj,ii) + I(nn)/(sqrt((rExt(ii)-rS(nn))^2 + (zExt(jj)-zS(nn))^2) + r0);
    end
   end
  end

%{
  %Add extra point sources if needed
   nSextr=3;
   Iextr(1)=-0.0007; zSextr(1)=0.4; rSextr(1)=0.2; r0extr(1)=0.1;
   Iextr(2)=0.0004; zSextr(2)=-0.6; rSextr(2)=-0.05; r0extr(2)=0.1;
   Iextr(3)=-0.0008; zSextr(3)=-0.15; rSextr(3)=-0.05; r0extr(3)=0.1;

   for nn = 1:nSextr
    for jj = 1:nhExt
     for ii = 1:nwExt
      psi_fit(jj,ii) = psi_fit(jj,ii) + Iextr(nn)/(sqrt((rExt(ii)-rSextr(nn))^2 + (zExt(jj)-zSextr(nn))^2) + r0extr(nn));
     end
    end
   end
%} 

  % Smooth the flux
  if s >= 0.
    [psi_smoothed psi_dct] = smoothn_mod(psi_fit, s);
  else
    [psi_smoothed psi_dct] = smoothn_mod(psi_fit);
  end

  % Write the smoothed flux
  fprintf(fd_out, '%16.9e%16.9e%16.9e%16.9e%16.9e\n', psi_smoothed');

  % Find the X point, given an initial guess (that might need to be adjusted)
  xpt = [0.75*rmin+0.25*rmax, 0.75*zmin+0.25*zmax];
  xpt = fminsearch(@(x) poloidal_Bmag(x, psi_dct', [rminExt, rmaxExt, zminExt, zmaxExt]), xpt);
  fprintf('The X point is at (R,Z) = (%f,%f)\n', xpt(1), xpt(2));

  % Write the dct coefficients
  fprintf(fd_dct, '%16.8f%16.8f%16.8f%16.8f%8d%8d\n', rminExt, rmaxExt, zminExt, zmaxExt, nwExt, nhExt);
  fprintf(fd_dct, '%16.8f%16.8f%16.8f%16.8f\n', rmaxis, zmaxis, xpt(1), xpt(2));
  fprintf(fd_dct, '%20.13e %20.13e %20.13e %20.13e %20.13e\n', psi_dct');


  % Make plots
  cmin = min(psi(:));
  cmax = max(psi(:));

  figure(1);
  h1 = contour(r,z,psi,80);
  caxis([cmin cmax]);
  axis image;

  figure(2);
  h2 = contour(rExt,zExt,psi_fit,80);
  caxis([cmin cmax]);
  axis image;

  cmin = min(psi_smoothed(:));
  cmax = max(psi_smoothed(:));

  figure(3);
  h4 = contour(rExt,zExt,psi_smoothed,80);
  caxis([cmin cmax]);
  axis image;
  
  figure(4);
  h3 = surf(r,z,(psi_fit - psi)/cmax);
  set(h3, 'EdgeColor', 'none');

  figure(5);
  h5 = surf(r,z,(psi_smoothed - psi)/cmax);
  set(h5, 'EdgeColor', 'none');


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

