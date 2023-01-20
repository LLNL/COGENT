%%%
% This script reads gfile (for geometry data), COGENT grid generated
% toroidal_flux_data (for q and Phi_norm(Psi_norm) verification) and 
% an inone file for plasma parameters. The script converts the original
% profiles as functions of rho=Sqrt(Phi_norm) to functions of Psi_norm, 
% and provides capabilites to extend across separatrix and to flatten at
% the inner core boundary.
%%%
clear all;

%Set filenames
inone_fname = 'inone.150142';
eqdsk_fname = 'g150142.02040';
cogent_tor_flux_filename = '../toroidal_flux_data';

% Read EQDSK file
fd_in  = fopen(eqdsk_fname,'r');

% Read and write the first line.  All we need is the problem size
tline = fgets(fd_in);
[header,count,errmsg,nextind] = sscanf(tline, '%s', 6);
n = sscanf(tline(nextind:end), '%d', 2);
nw = n(1);
nh = n(2);

% Read the second and third lines containing the geometry data
% we want
tline = fgets(fd_in);
geom1 = sscanf(tline, '%g', 5);
tline = fgets(fd_in);
geom2 = sscanf(tline, '%g', 5);

% Read and write some lines we don't care about
for i=1:2
  tline = fgets(fd_in);
end

% Read RBtor profile
RBtor_profile = fscanf(fd_in, '%g', nw);

% Read pressure and GS-related parameters
pres = fscanf(fd_in, '%g', nw);
ffprim = fscanf(fd_in, '%g', nw);
pprime = fscanf(fd_in, '%g', nw);

% Read the flux
psi_orig = fscanf(fd_in, '%g', [nw nh])';

% Read the safety factor
safety_factor = fscanf(fd_in, '%g', nw);

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

% Get some field data
psi_axis = geom2(3);
psi_sep = geom2(4);
Bvac_centr = geom2(5);
RBtor_vac = Bvac_centr*rcentr;

% Evaluate toroidal flux as function of a poloidal flux by q-integration
psi_nodal_grid = linspace(psi_axis, psi_sep, nw);
psi_norm_nodal_grid = (psi_nodal_grid(:)-psi_axis)/(psi_sep-psi_axis);
toroidal_flux_upper_face = zeros(1,nw-1);
psi_norm_upper_face = zeros(1,nw-1);
for n=1:nw-1
    safety_factor_cc = 0.5*(safety_factor(n)+safety_factor(n+1));
    dpsi = psi_nodal_grid(n+1)-psi_nodal_grid(n);
    psi_norm_upper_face(n)=psi_norm_nodal_grid(n+1);
    if n==1
        toroidal_flux_upper_face(n) = 2*pi*safety_factor_cc * dpsi;
    else 
        toroidal_flux_upper_face(n) = toroidal_flux_upper_face(n-1)+2*pi*safety_factor_cc * dpsi;
    end
end

norm_toroidal_flux_upper_face = toroidal_flux_upper_face(:)/toroidal_flux_upper_face(nw-1);

%%% Read toroidal flux and safety factor data from COGENT grid generator
Tcogent=readtable(cogent_tor_flux_filename);

%%% Plot Phi_norm and safety factor for eqdsk and cogent grid generator calcualtions
%%% N.B. for the grid calculation safety factor is cell-centered, but
%%% poloidal-flux is upper-faced. So, formally there is first order dr
%%% offset, but given a really small cell size, we don't care about it. 
f1 = figure(1);
hold on
yyaxis left
ylabel('Normalized toroidal flux, \Phi_n = \Phi/\Phi_{sep}')
plot(psi_norm_upper_face,norm_toroidal_flux_upper_face);
plot(Tcogent.(2)(:),Tcogent.(1)(:));
yyaxis right
ylabel('safety factor, q')
plot(linspace(0,1,nw), safety_factor);
plot(Tcogent.(2)(:),abs(RBtor_vac)*Tcogent.(3)(:));
hold off;
title('COGENT vs eqdsk')
xlabel('normalized poloidal flux')
legend({'\Phi_n eqdsk', '\Phi_n COGENT', 'q COGENT', 'q eqdsk'},'Location','northwest')
ax = gca;
ax.FontSize = 13;


%%% Plot RBtor profile
f3 = figure(3);
plot(linspace(0,1,nw),RBtor_profile);
yline(RBtor_vac, '-', 'RBtor vaccum');
title('RBtor profile (T \cdot m)')
xlabel('normalized poloidal flux')
ax = gca;
ax.FontSize = 13;

%%% Read inone file for plasma parameters
[rho_input, ne_input, Ti_input, Te_input, Zeff_input] = readInOneFile(inone_fname);

%%% Plot onetwo plasma profiles as function of rho
f2 = figure(2);
hold on
plot(rho_input,ne_input/10^13);
plot(rho_input,Ti_input);
plot(rho_input, Te_input);
hold off;

title('Core profiles')
xlabel('\rho = (\Phi/\Phi_{sep})^{1/2}')
ylabel('Normalized data')
legend({'Ne (10^{13}/cm^3)', 'Ti (keV)', 'Te (keV)'},'Location','southwest')
ax = gca;
ax.FontSize = 13;

%%% Post process inone data
%%% phi and psi fluxes are assumed normalized here

psi_core = 0.85; % valid inner core bndry
psi_sol = 1.05;  % valid outer sol bndry

psi_core_ghost = psi_core - 0.05; % inner core ghost location
psi_sol_ghost = psi_sol + 0.05;   %% outer sol ghost location

psi_match_core = 0.89; % matching location near core bndry
psi_match_sol = 1.0;   % matching location near separatrix

Npts = 60; % number of point for the output psi grid
psi_output = linspace(psi_core_ghost, psi_sol_ghost, Npts);

% Convert rho to psi_norm
phi_input(:) = rho_input(:).^2; 
psi_input = interp1(norm_toroidal_flux_upper_face, psi_norm_upper_face, phi_input, 'spline');

% Setting piecewise polinomial spline for density
density = extrapolateProfile(psi_input, ne_input, psi_output, psi_core, psi_sol, psi_match_core, psi_match_sol, 1.0, true);
ion_temperature = extrapolateProfile(psi_input, Ti_input, psi_output, psi_core, psi_sol, psi_match_core, psi_match_sol, 1.0, true);
electron_temperature = extrapolateProfile(psi_input, Te_input, psi_output, psi_core, psi_sol, psi_match_core, psi_match_sol, 4.0, false);
Zeff = extrapolateProfile(psi_input, Zeff_input, psi_output, psi_core, psi_sol, psi_match_core, psi_match_sol, 4.0, false);

% Getting normalization
pp_ne=spline(psi_input, ne_input);
ne_norm = ppval(pp_ne, psi_core);

pp_Ti=spline(psi_input, Ti_input);
Ti_norm = ppval(pp_Ti, psi_core);

%%% Plot plasma profiles %%%%
f4 = figure(4);
psi_core_grid = linspace(psi_core_ghost,1.0,50);
hold on
plot(psi_core_grid,ppval(pp_ne,psi_core_grid)/ne_norm, '*');
plot(psi_output,density/ne_norm);
plot(psi_core_grid,ppval(pp_Ti,psi_core_grid)/Ti_norm, '*');
plot(psi_output,ion_temperature/Ti_norm);
plot(psi_output,electron_temperature/Ti_norm);
hold off;

title('Edge profiles with extrapolation')
xlabel('normalized poloidal flux')
ylabel('Normalized data')

legend({'Ne orig','Ne pp', 'Ti orig', 'Ti pp', 'Te pp'},'Location','southwest')
ax = gca;
ax.FontSize = 13;

txt1 = append('T_0(keV) = ',num2str(Ti_norm,3)) ;
txt2 = append('N_0(10^{13}/cm^3) = ',num2str(ne_norm/10^13,3)) ;
text(1,1,{txt1,txt2},'FontSize',14)

%%% Comput parallel conductivity 
pp_q = spline(linspace(0,1,nw), safety_factor);
pp_Zeff = spline(psi_output, Zeff);
for i=1:Npts
    x=psi_output(i);
    q(i) = ppval(pp_q, x);
    if (x>0.85) 
        q(i) = ppval(pp_q, 0.85); % limit q by psi=0.95 value
        Zeff(i) = ppval(pp_Zeff, 0.85); % limit Zeff by psi=0.95 value
    end
end

[sigma, psi_crit] = computeParallelConductivity(psi_output, density, electron_temperature, Zeff, q, rcentr, ne_norm, Ti_norm);

%%% Print plasma profile data
ne(:,1)=psi_output;
ne(:,2)=density/ne_norm;
writematrix(ne,'ne.txt', 'Delimiter','tab') 

Ti(:,1)=psi_output;
Ti(:,2)=ion_temperature/Ti_norm;
writematrix(Ti,'Ti.txt', 'Delimiter','tab') 

Te(:,1)=psi_output;
Te(:,2)=electron_temperature/Ti_norm;
writematrix(Te,'Te.txt', 'Delimiter','tab') 

sig(:,1)=psi_output;
sig(:,2)=sigma;
writematrix(sig,'sigma.txt', 'Delimiter','tab') 

%%% Internal functions
function profile = extrapolateProfile(psi_input, data, psi_output, psi_core, psi_sol, psi_match_core, psi_match_sol, sol_smoothing_factor, core_smoothing)


  delta_core = (psi_core - psi_match_core)/1.0;
  delta_sol = (psi_sol - psi_match_sol)/sol_smoothing_factor;

  % Setting piecewise polinomial (pp) spline 
  ppf=spline(psi_input, data);
  
  % Getting pp derivative 
  ppdf = ppf;
  ppdf.order=ppdf.order-1;
  ppdf.coefs=ppdf.coefs(:,1:end-1).*(ppdf.order:-1:1);

  % Construct f0+f1*tanch((-x+xmatch)/delta)
  f0_core = ppval(ppf, psi_match_core);
  f1_core = -ppval(ppdf,psi_match_core)*delta_core;

  f0_sol = ppval(ppf, psi_match_sol);
  f1_sol = -ppval(ppdf, psi_match_sol)*delta_sol;

  % Created smoothed and extrapolated function
  Npt = length(psi_output);
  profile = zeros(1,Npt);
  for i=1:Npt
    x=psi_output(i);
    if x <= psi_match_core 
        if (core_smoothing == true) 
            profile(i) = f0_core + f1_core * tanh((-x+psi_match_core)/delta_core); 
        else 
            profile(i) = ppval(ppf, x);
        end
    elseif (x > psi_match_core) && (x <= psi_match_sol)
        profile(i) = ppval(ppf, x);
    else
        profile(i) = f0_sol + f1_sol * tanh((-x+psi_match_sol)/delta_sol); 
    end
  end
end

function [sigma, psi_crit] = computeParallelConductivity(psi_output, ne, Te, Zeff, q, R0, n_norm, T_norm)
  
    %%% Computes normalized sigma and  
    %%% the value of psi at which nu_ei * qR0/VTe = 1
    %%% N.B. q profile that is passed here is flattend at some core value

    ech = 4.8032*10^(-10);
    mp = 1.6726*10^(-24);
    me = 9.1*10^(-28);

    Tcgs_norm = T_norm*1000*1.602*10^(-12);
    Lcgs_norm = 100; %1m
    Vcgs_norm = sqrt(Tcgs_norm/mp); 

    R0_cgs = R0*100;
    
    % Set psi_crit to the left boundary
    psi_crit = psi_output(1);

    Npt = length(psi_output);
    for i=1:Npt
        TeeV = Te(i) * 1000;
        LnLambdaEl = 23.5 - log(TeeV^(-5/4)*ne(i)^(1/2))-sqrt(10^(-5) + 1/16*(log(TeeV) - 2)^2);
        tauBragE = 3.44*10^(5)*TeeV^(3/2)/ne(i)/Zeff(i)/LnLambdaEl;

        Vte = sqrt(TeeV*1.602*10^(-12)/me);
        if (tauBragE>q(i)*R0_cgs/Vte)
            tauBragE = q(i)*R0_cgs/Vte;
            psi_crit = psi_output(i);
        end    

        sigma_unnorm = ech^2 * ne(i) *tauBragE/(0.51 * me);
        sigma(i) = sigma_unnorm * Tcgs_norm / (ech^2 *Vcgs_norm * n_norm * Lcgs_norm);
    end

end

function [RENEIN, ENEIN , TIIN, TEIN, ZEFFIN] = readInOneFile(inone_fname)
    fileID = fopen(inone_fname);

    %%%  Read in each element of the file and store as cell array
    C = textscan(fileID,'%s');
    NC = length(C{1});
 
    %%%  Find the number of data points for the 1D arrays
    NJENE_str = 'NJENE';

    for i=1:NC
        this_str = C{1}{i};
        if(strcmp(NJENE_str,this_str))
            index_NJENE = i+2;
            break;
        end
    end
 
    Ndata = str2num(C{1}{index_NJENE});

    %%% Get radial coordinates (normalized toroidal flux)
    data_str = 'RENEIN(1,1)';

    for i=index_NJENE:NC
        this_str = C{1}{i};
        if(strcmp(data_str,this_str))
            index_data = i+2;
            break;
        end
    end

    RENEIN = zeros(1,Ndata);
    for n=1:Ndata
        i = index_data + n - 1;
        RENEIN(n) = str2num(C{1}{i});
    end

    %%% Getn ion density m^-3
    data_str = 'ENEIN(1,1)';

    for i=index_data:NC
        this_str = C{1}{i};
        if(strcmp(data_str,this_str))
            index_data = i+2;
            break;
        end
    end

    ENEIN = zeros(1,Ndata);
    for n=1:Ndata
        i = index_data + n - 1;
        ENEIN(n) = str2num(C{1}{i});
    end

    %%% Getn electron temperature, keV
    data_str = 'TEIN(1,1)';

    for i=index_data:NC
        this_str = C{1}{i};
        if(strcmp(data_str,this_str))
            index_data = i+2;
            break;
        end
    end

    TEIN = zeros(1,Ndata);
    for n=1:Ndata
        i = index_data + n - 1;
        TEIN(n) = str2num(C{1}{i});
    end

    %%% Getn ion temperature, keV
    data_str = 'TIIN(1,1)';

    for i=index_data:NC
        this_str = C{1}{i};
        if(strcmp(data_str,this_str))
            index_data = i+2;
            break;
        end
    end

    TIIN = zeros(1,Ndata);
    for n=1:Ndata
        i = index_data + n - 1;
        TIIN(n) = str2num(C{1}{i});
    end

    %%% Getn Zeff
    data_str = 'ZEFFIN(1,1)';

    for i=index_data:NC
        this_str = C{1}{i};
        if(strcmp(data_str,this_str))
            index_data = i+2;
            break;
        end
    end

    ZEFFIN = zeros(1,Ndata);
    for n=1:Ndata
        i = index_data + n - 1;
        ZEFFIN(n) = str2num(C{1}{i});
    end

    %%% Get logarithmic gradients
    for n=1:Ndata-1
        KappaT(n) = 2*(TIIN(n+1)-TIIN(n))/((TIIN(n+1)+TIIN(n)));
        KappaN(n) = 2*(ENEIN(n+1)-ENEIN(n))/((ENEIN(n+1)+ENEIN(n)));
    end
end
