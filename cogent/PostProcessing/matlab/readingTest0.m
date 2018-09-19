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

thisPath1 = './D3D_multiBlock/plt_density_plots/';
thisFile1 = [thisPath1,'plt.1.hydrogen.density0000.2d.hdf5'];
%
%thisPath2 = './D3D_multiBlock/plt_efield_plots/';
%thisFile2 = [thisPath2,'plt.efield0000.2d.hdf5'];
%
thisPath2 = './D3D_multiBlock/plt_potential_plots/';
thisFile2 = [thisPath2,'plt.potential0390.2d.hdf5'];
%%%   load the data
%
%data1 = import2Ddata_singleFile(thisFile1);
data2 = import2Ddata_singleFile(thisFile2);


%%%
%
fileinfo1 = hdf5info(thisFile1);
GH1 = fileinfo1.GroupHierarchy;
ghost1 = h5readatt(thisFile1,'/level_0/data_attributes','ghost');
outputGhost1 = h5readatt(thisFile1,'/level_0/data_attributes','outputGhost');
nComps1 = h5readatt(thisFile1,'/','num_components');
prob_domain1 = h5readatt(thisFile1,'/level_0','prob_domain');
%ref_ratio1 = h5readatt(thisFile1,'/level_0','ref_ratio');
vecData1  = hdf5read(thisFile1,'/level_0/data:datatype=0');
offsets1  = hdf5read(thisFile1,'/level_0/data:offsets=0');
boxes1    = hdf5read(thisFile1,'/level_0/boxes');
size(vecData1);
%
thisFileMap1 = [thisFile1(1:end-5),'.map.hdf5'];
fileinfoMap1 = hdf5info(thisFileMap1);
vecMap1      = hdf5read(thisFileMap1,'/level_0/data:datatype=0');
offsetsMap1  = hdf5read(thisFileMap1,'/level_0/data:offsets=0');
ghostMap1 = h5readatt(thisFileMap1,'/level_0/data_attributes','ghost');
boxesMap1    = hdf5read(thisFileMap1,'/level_0/boxes');

size(vecMap1);


%%%
%
fileinfo2 = hdf5info(thisFile2);
GH2 = fileinfo2.GroupHierarchy;
ghost2 = h5readatt(thisFile2,'/level_0/data_attributes','ghost');
outputGhost2 = h5readatt(thisFile2,'/level_0/data_attributes','outputGhost');
nComps2 = h5readatt(thisFile2,'/','num_components');
prob_domain2 = h5readatt(thisFile2,'/level_0','prob_domain');
%ref_ratio2 = h5readatt(thisFile2,'/level_0','ref_ratio');
vecData2  = hdf5read(thisFile2,'/level_0/data:datatype=0');
offsets2  = hdf5read(thisFile2,'/level_0/data:offsets=0');
boxes2    = hdf5read(thisFile2,'/level_0/boxes');
size(vecData2);
%
thisFileMap2 = [thisFile2(1:end-5),'.map.hdf5'];
fileinfoMap2 = hdf5info(thisFileMap2);
vecMap2      = hdf5read(thisFileMap2,'/level_0/data:datatype=0');
offsetsMap2  = hdf5read(thisFileMap2,'/level_0/data:offsets=0');
ghostMap2 = h5readatt(thisFileMap2,'/level_0/data_attributes','ghost');
boxesMap2    = hdf5read(thisFileMap2,'/level_0/boxes');
size(vecMap2);





