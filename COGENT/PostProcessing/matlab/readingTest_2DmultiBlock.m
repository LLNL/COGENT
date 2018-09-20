%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   reading test for 2D data with and without multi-block
%%%
%%%   There are 2 different 2D reader functions
%%%   1) import2Ddata_singleFile(thisFile), which takes a single file name
%%%      as an input and returns the data from that single file
%%%   2) import2Ddata(thisPath,thisSpecies), which takes the path to a folder
%%%      containing any number of outputs of the same type corresponding
%%%      to different time outputs and loads them all into a single
%%%      structure. The second argument is the species number. This
%%%      arguement defaults to 1 if not specified.
%%%
%%%   NOTE: The folder with the data is assumed to also contain the .map files
%%%
%%%%%%%%%%%%%%%%%%%%%


clear all;
addpath('~angus1/Programs/COGENT_matlabTools/');


%%%   specify folder path and file name (choose one of the three sets
%%%   below)
%
%thisPath = './Torus_singleBlock/plt_vpar_mu_plots/';
%thisFile = [thisPath,'plt.1.hydrogenA.vpar_mu0000.2d.hdf5']; species = 1;
%
%thisPath = './Torus_singleBlock/plt_density_plots/';
%thisFile = [thisPath,'plt.1.hydrogenA.density0000.2d.hdf5']; species = 1;
%
thisPath = './D3D_multiBlock/plt_density_plots/';
thisFile = [thisPath,'plt.1.hydrogen.density0000.2d.hdf5']; species = 1;
%
thisPath = './D3D_multiBlock/plt_efield_plots/';
thisFile = [thisPath,'plt.efield0390.2d.hdf5']; species = 0;
%
thisPath = './D3D_multiBlock/plt_potential_plots/';
thisFile = [thisPath,'plt.potential0390.2d.hdf5']; species = 0;

%%%   load the data
%
data = import2Ddata_singleFile(thisFile);
data2 = import2Ddata(thisPath,species);


nGX = data.numGhostX;
nGY = data.numGhostY;

%%%   plot one of the grid variables and the data on unmapped output
%%%   Note that multiblock output has a lot of zeros in between blocks.
%
close(figure(1));
f1 = figure(1); set(f1,'position',[70 1 1100 800]);
%
subplot(1,2,1);
pcolor(data.Xce); title('unmapped X'); colorbar; axis('square');
%
subplot(1,2,2);
pcolor(data.Yce); title('unmapped Y'); colorbar; axis('square');


%%%   plot data on unmapped and mapped grid
%
close(figure(2));
f2 = figure(2); set(f2,'position',[70 1 1100 800]);
%
subplot(1,2,1);
pcolor(data.Fcc(:,:,1)); title('unmapped data'); colorbar; axis('square');
%
numBlocks = data.numBlocks;
for b=1:numBlocks
    Xce = data.block(b).Xce;
    Yce = data.block(b).Yce;
    Fce = data.block(b).Fce;
   % Fce = squeeze(data2.block(b).Fce(:,:,:,1));
    figure(2); hold on;
    subplot(1,2,2); 
    hold on; pcolor(Xce(nGX+1:end-nGX,nGY+1:end-nGY), ...
                    Yce(nGX+1:end-nGX,nGY+1:end-nGY), ...
                    Fce(nGX+1:end-nGX,nGY+1:end-nGY,1)); box on;
   % end
    box on; axis('equal'); colorbar;
    title('mapped data');
end


fileinfo = hdf5info(thisFile);
GH = fileinfo.GroupHierarchy;
vecData  = hdf5read(thisFile,'/level_0/data:datatype=0');
ghost = h5readatt(thisFile,'/level_0/data_attributes','ghost');
nComps = h5readatt(thisFile,'/','num_components');
size(vecData);
%
thisFileMap = [thisFile(1:end-5),'.map.hdf5'];
fileinfoMap = hdf5info(thisFileMap);
vecMap      = hdf5read(thisFileMap,'/level_0/data:datatype=0');
ghostMap = h5readatt(thisFileMap,'/level_0/data_attributes','ghost');
size(vecMap);
%figure(11); plot(vecData);





