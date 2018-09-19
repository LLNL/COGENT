%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   reading test for 4D data (x,z,vpar,mu) with and without multi-block
%%%
%%%   There are 2 different 4D reader functions
%%%   1) import4Ddata_singleFile(thisFile), which takes a single file name
%%%      as an input and returns the data from that single file
%%%   2) import4Ddata(thisPath,thisSpecies), which takes the path to a folder
%%%      containing any number of outputs of the same type corresponding
%%%      to different time outputs and loads them all into a single
%%%      structure. The second argument is the species number. This
%%%      arguement defaults to 1 if not specified.
%%%
%%%   NOTE: The folder with the data is assumed to also contain the .map files
%%%
%%%   NOTE: It is assumed multiblock on for real space dims
%%%
%%%%%%%%%%%%%%%%%%%%%
clear all;
addpath('~angus1/Programs/COGENT_matlabTools/');


%%%   specify folder path (choose one of the two sets below)
%
thisPath = './Torus_singleBlock/plt_dfn_plots/';
thisFile = [thisPath,'plt.1.hydrogenA.dfn0000.4d.hdf5'];
data = import4Ddata(thisPath);
data2 = import4Ddata_singleFile(thisFile);
%
%thisPath = './D3D_multiBlock/plt_dfn_plots/';
%thisFile = [thisPath,'plt.1.hydrogen.dfn0000.4d.hdf5'];
%data = import4Ddata(thisPath);
%data2 = import4Ddata_singleFile(thisFile);


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
vparIndex = round(length(data.block(1).Xce(1,1,:,1))/2);
muIndex = 1; %round(length(data.block(1).Xce(1,1,1,:))/2);
%
subplot(2,2,1);
numBlocks = data.numBlocks;
for b=1:numBlocks
%hold on; pcolor(data.block(b).Xce(:,:,vparIndex,muIndex),data.block(b).X2ce(:,:,vparIndex,muIndex),data.block(b).Fce(:,:,vparIndex,muIndex,1)); 
hold on; pcolor(data2.block(b).Xce(:,:,vparIndex,muIndex),data2.block(b).X2ce(:,:,vparIndex,muIndex),data2.block(b).Fce(:,:,vparIndex,muIndex)); 
end
box on;
title(['data at vpar=',num2str(data.X3ce(1,1,vparIndex,muIndex)),', \mu=',num2str(data.X4ce(1,1,vparIndex,muIndex))]); 
colorbar; axis('equal'); xlabel('r'); ylabel('z');
%
subplot(2,2,2);
pcolor(squeeze(data.block(1).X3ce(1,1,:,:)),squeeze(data.block(1).X4ce(1,1,:,:)), ...
       squeeze(data.block(1).Fce(1,1,:,:,1))); 
title(['data at r=',num2str(data.block(1).Xce(1,1,1,1)),', z=',num2str(data.block(1).X2ce(1,1,1,1))]); 
colorbar; axis('square'); xlabel('vpar'); ylabel('\mu');
%
subplot(2,2,3);
pcolor(squeeze(data.block(1).Xce(:,1,vparIndex,:)),squeeze(data.block(1).X4ce(:,1,vparIndex,:)), ...
       squeeze(data.block(1).Fce(:,1,vparIndex,:,1))); 
title(['data at \theta=',num2str(data.block(1).X2ce(1,1,1,1)),', vpar=',num2str(data.block(1).X3ce(1,1,vparIndex,1))]); 
colorbar; axis('square'); xlabel('r'); ylabel('\mu');
%
subplot(2,2,4);
pcolor(squeeze(data.block(1).X3ce(1,:,:,muIndex)),squeeze(data.block(1).X2ce(1,:,:,muIndex)), ...
       squeeze(data.block(1).Fce(1,:,:,muIndex,1))); 
title(['data at r=',num2str(data.block(1).X2ce(1,1,1,1)),', \mu=',num2str(data.block(1).X4ce(1,1,1,muIndex))]); 
colorbar; axis('square'); xlabel('vpar'); ylabel('z');


