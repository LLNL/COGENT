%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   reading single file in 4D (x,y,vpar,mu)
%%%
%%%%%%%%%%%%%%%%%%%%%
clear all;
addpath('~angus1/Programs/COGENT_matlabTools/');


thisPath = './Torus_singleBlock/plt_dfn_plots/';
data = import4Ddata(thisPath);

% thisPath = './D3D_multiBlock/plt_dfn_plots/';
% data = import4Ddata(thisPath);

%%%   plot one of the grid variables and the data on unmapped output
%%%   Note that multiblock output has a lot of zeros in between blocks.
%
close(figure(1));
f1 = figure(1); set(f1,'position',[1220 100 1250 1150]);
%
subplot(2,2,1);
pcolor(data.Xce(:,:,1,1)'); title('unmapped X'); colorbar; axis('square');
xlabel('radial index'); ylabel('poloidal index');
%
subplot(2,2,2);
pcolor(data.X2ce(:,:,1,1)); title('unmapped Y'); colorbar; axis('square');
xlabel('poloidal index'); ylabel('radial index');
%
subplot(2,2,3);
pcolor(squeeze(data.X3ce(:,1,:,1))); title('unmapped vpar'); colorbar; axis('square');
xlabel('vpar index'); ylabel('radial index');
%
subplot(2,2,4);
pcolor(squeeze(data.X4ce(:,1,1,:))); title('unmapped mu'); colorbar; axis('square');
xlabel('mu index'); ylabel('radial index');



%%%   plot several 2D phase space contours of 4D data
%
close(figure(2));
f2 = figure(2); set(f2,'position',[1220 100 1250 1150]);
vparIndex = round(length(data.Fce(1,1,:,1))/2);
muIndex = 1; %round(length(data.Fce(1,1,1,:))/2);
%
subplot(2,2,1);
pcolor(data.Xce(:,:,vparIndex,muIndex),data.X2ce(:,:,vparIndex,muIndex),data.Fce(:,:,vparIndex,muIndex,1)); 
title(['data at vpar=',num2str(data.X3ce(1,1,vparIndex,muIndex)),', \mu=',num2str(data.X4ce(1,1,vparIndex,muIndex))]); 
colorbar; axis('square'); xlabel('r'); ylabel('z');
%
subplot(2,2,2);
pcolor(squeeze(data.X3ce(1,1,:,:)),squeeze(data.X4ce(1,1,:,:)), ...
       squeeze(data.Fce(1,1,:,:,1))); 
title(['data at r=',num2str(data.Xce(1,1,1,1)),', z=',num2str(data.X2ce(1,1,1,1))]); 
colorbar; axis('square'); xlabel('vpar'); ylabel('\mu');
%
subplot(2,2,3);
pcolor(squeeze(data.Xce(:,1,vparIndex,:)),squeeze(data.X4ce(:,1,vparIndex,:)), ...
       squeeze(data.Fce(:,1,vparIndex,:,1))); 
title(['data at \theta=',num2str(data.X2ce(1,1,1,1)),', vpar=',num2str(data.X3ce(1,1,vparIndex,1))]); 
colorbar; axis('square'); xlabel('r'); ylabel('\mu');
%
subplot(2,2,4);
pcolor(squeeze(data.X3ce(1,:,:,muIndex)),squeeze(data.X2ce(1,:,:,muIndex)), ...
       squeeze(data.Fce(1,:,:,muIndex,1))); 
title(['data at r=',num2str(data.X2ce(1,1,1,1)),', \mu=',num2str(data.X4ce(1,1,1,muIndex))]); 
colorbar; axis('square'); xlabel('vpar'); ylabel('z');


