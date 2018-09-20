function smoothndemo

%SMOOTHNDEMO Demonstration for the use of SMOOTHN.
%   SMOOTHNDEMO provides four examples that illustrate data smoothing with
%   SMOOTHN. The figures correspond to the Figures 2 to 5 of the paper
%   entitled "Fast and robust smoothing of gridded data in one and higher
%   dimensions with occurrence of missing values", written by Damien Garcia
%   and published in Computational Statistics and Data Analysis (2009).
%

%% Figure 2
% Automatic smoothing of 2-D data with missing values
% ---
figure(2)
% Original data
n = 300;
y0 = peaks(n);
miny0 = min(y0(:)); maxy0 = max(y0(:));
% Add Gaussian noise
yn = y0 + randn(size(y0));
subplot(221), imagesc(yn), title('Noisy data')
caxis([miny0 maxy0]), axis equal off
% Lose 1/2 of data
I = randperm(n^2);
yn(I(1:n^2*1/2)) = NaN;
% Create a hole
yn(100:155,75:125) = NaN;
subplot(222), imagesc(yn), title('Noisy data with missing values')
caxis([miny0 maxy0]), axis equal off
% Smooth data
z = smoothn(yn,'MaxIter',50);
subplot(223), imagesc(z), title('Smoothed data')
caxis([miny0 maxy0]), axis equal off
% Absolute errors
subplot(224), imagesc(z-y0), title('Absolute errors')
caxis([miny0 maxy0]), axis equal off

%% Figure 3 (see also example 1 in Appendix II)
% Non robust versus robust smoothing (1-D)
% ---
figure(3)
% Original data + white noise
x = linspace(0,100,2^8);
y = cos(x/10)+(x/50).^2 + randn(size(x))/10;
% Add outliers
y([70 75 80]) = [5.5 5 6];
% Non robust smoothing
z = smoothn(y);
subplot(121), plot(x,y,'r.',x,z,'k','LineWidth',2)
axis square, title('Non robust')
% Robust smoothing
z = smoothn(y,'robust');
subplot(122), plot(x,y,'r.',x,z,'k','LineWidth',2)
axis square, title('Robust')


%% Figure 4
% Automatic smoothing of 3-D data
% ---
figure(4), colormap(hsv)
[x,y,z] = meshgrid(-2:.2:2,-2:.2:2,-2:.2:2);
% Visualization specs
az = -24; el = 23;
xslice = [-0.8,1]; yslice = 2; zslice = [-2,0];
% Original data
v = x.*exp(-x.^2-y.^2-z.^2);
% Noisy data
vn = v + randn(size(v))*0.06;
subplot(121)
slice(x,y,z,vn,xslice,yslice,zslice,'cubic')
title('Noisy data')
caxis([min(v(:)) max(v(:))]), view(az,el), axis off equal
% Smoothed data
vs = smoothn(vn);
subplot(122)
slice(x,y,z,vs,xslice,yslice,zslice,'cubic')
title('Smoothed data')
caxis([min(v(:)) max(v(:))]), view(az,el), axis off equal


%% Figure 5
% Automatic smoothing of a noisy cardioid
% ---
figure(5)
% Cardioid
t = linspace(0,2*pi,1000);
x = 2*cos(t).*(1-cos(t));
y = 2*sin(t).*(1-cos(t));
% Noisy cardioid
xn = x + randn(size(x))*0.15;
yn = y + randn(size(y))*0.15;
% Smoothed cardioid
z = smoothn(complex(xn,yn));
plot(xn,yn,'r.',real(z),imag(z),'k','linewidth',2)
axis equal tight
title('Smoothing of a noisy cardioid')

