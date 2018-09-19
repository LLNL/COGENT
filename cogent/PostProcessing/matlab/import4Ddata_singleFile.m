function [DataStr] = import4Ddata_singleFile(fileName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   import 4D data from single COGENT output file
%%%
%%%   input parameters:
%%%   fileName = 'string/to/plt_something_plots/thisFile' 
%%%
%%%   output parameter is a structure with:
%%%   block structure with Xce, X2ce, X3ce, X4ce, and Fce for each block
%%%   time found in output file
%%%   cell-edge grid (Xce and X2ce, X3ce, X4ce) for all blocks
%%%   cell-center data (Fcc) for all blocks
%%%   number of processors
%%%   number of blocks
%%%   number of ghost cells at each end in X
%%%   number of ghost cells at each end in Y
%%%
%%%   block structure has local Xce, X2ce, X3ce, X4ce, and Fce for each block, where 
%%%   Fce is cell centered data extened by one in each direction
%%%   to be compatible with Xce and Yce for pcolor and contour plots
%%%
%%%   Note that the function value is defined at cell-center
%%%   while the grid is defined at cell-edge
%%%
%%%   March 28, 2018
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thisFile = fileName;
fileinfo = hdf5info(thisFile);
time = h5readatt(thisFile,'/level_0','time');
%  dx = h5readatt(thisFile,'/level_0','dx');
prob_domain = h5readatt(thisFile,'/level_0','prob_domain');
ghost = h5readatt(thisFile,'/level_0/data_attributes','ghost');
outputGhost = h5readatt(thisFile,'/level_0/data_attributes','outputGhost');
nx = double(prob_domain.hi_i-prob_domain.lo_i+1);
nx2 = double(prob_domain.hi_j-prob_domain.lo_j+1);
nx3 = double(prob_domain.hi_k-prob_domain.lo_k+1);
nx4 = double(prob_domain.hi_l-prob_domain.lo_l+1);
procs = hdf5read(thisFile,'/level_0/Processors');
numProcs = length(procs);


%%%   loop over time and load the data
%
data0out = zeros(nx,nx2,nx3,nx4);
vecData  = hdf5read(thisFile,'/level_0/data:datatype=0');
offsets  = hdf5read(thisFile,'/level_0/data:offsets=0');
boxes    = hdf5read(thisFile,'/level_0/boxes');
for iP = 1:numProcs
    lo_i(iP) = boxes(iP).Data{1}+1-ghost.intvecti;
    lo_j(iP) = boxes(iP).Data{2}+1-ghost.intvectj;
    lo_k(iP) = boxes(iP).Data{3}+1-ghost.intvectk;
    lo_l(iP) = boxes(iP).Data{4}+1-ghost.intvectl;
    hi_i(iP) = boxes(iP).Data{5}+1+ghost.intvecti;
    hi_j(iP) = boxes(iP).Data{6}+1+ghost.intvectj;
    hi_k(iP) = boxes(iP).Data{7}+1+ghost.intvectk;
    hi_l(iP) = boxes(iP).Data{8}+1+ghost.intvectl;
end
    
%%%   some indices are negative, so shift for MATLAB
%
min_lo_i = min(lo_i);
lo_i = lo_i-double(min_lo_i)+1;
hi_i = hi_i-double(min_lo_i)+1;
min_lo_j = min(lo_j);
lo_j = lo_j-double(min_lo_j)+1;
hi_j = hi_j-double(min_lo_j)+1;
min_lo_k = min(lo_k);
lo_k = lo_k-double(min_lo_k)+1;
hi_k = hi_k-double(min_lo_k)+1;
min_lo_l = min(lo_l);
lo_l = lo_l-double(min_lo_l)+1;
hi_l = hi_l-double(min_lo_l)+1;    
    
    
thisFileMap = [fileName(1:end-5),'.map.hdf5'];
fileinfoMap = hdf5info(thisFileMap);
vecMap      = hdf5read(thisFileMap,'/level_0/data:datatype=0');
offsetsMap  = hdf5read(thisFileMap,'/level_0/data:offsets=0');
boxesMap    = hdf5read(thisFileMap,'/level_0/boxes');
lo_i_Map = zeros(1,numProcs);
lo_j_Map = zeros(1,numProcs);
lo_k_Map = zeros(1,numProcs);
lo_l_Map = zeros(1,numProcs);
hi_i_Map = zeros(1,numProcs);
hi_j_Map = zeros(1,numProcs);
hi_k_Map = zeros(1,numProcs);
hi_l_Map = zeros(1,numProcs);
for iP = 1:numProcs
    lo_i_Map(iP) = boxesMap(iP).Data{1}+1-ghost.intvecti;
    lo_j_Map(iP) = boxesMap(iP).Data{2}+1-ghost.intvectj;       
    lo_k_Map(iP) = boxesMap(iP).Data{3}+1-ghost.intvectk;
    lo_l_Map(iP) = boxesMap(iP).Data{4}+1+ghost.intvectl; 
    hi_i_Map(iP) = boxesMap(iP).Data{5}+1+ghost.intvecti;
    hi_j_Map(iP) = boxesMap(iP).Data{6}+1+ghost.intvectj;
    hi_k_Map(iP) = boxesMap(iP).Data{7}+1+ghost.intvectk;
    hi_l_Map(iP) = boxesMap(iP).Data{8}+1+ghost.intvectl;
end
min_lo_i_Map = min(lo_i_Map);
lo_i_Map = lo_i_Map-double(min_lo_i_Map)+1;
hi_i_Map = hi_i_Map-double(min_lo_i_Map)+1;
min_lo_j_Map = min(lo_j_Map);
lo_j_Map = lo_j_Map-double(min_lo_j_Map)+1;
hi_j_Map = hi_j_Map-double(min_lo_j_Map)+1;
min_lo_k_Map = min(lo_k_Map);
lo_k_Map = lo_k_Map-double(min_lo_k_Map)+1;
hi_k_Map = hi_k_Map-double(min_lo_k_Map)+1;
min_lo_l_Map = min(lo_l_Map);
lo_l_Map = lo_l_Map-double(min_lo_l_Map)+1;
hi_l_Map = hi_l_Map-double(min_lo_l_Map)+1;

    
%%%   Note that boxes and boxesMap are always the same, while   %%%
%%%   data is at cell center and grid at cell-edges             %%%
    

%%%   map to reshaped grid
%
X = zeros(nx+1,nx2+1,nx3+1,nx4+1);   % at cell edge
X2 = zeros(nx+1,nx2+1,nx3+1,nx4+1);   % at cell edge
X3 = zeros(nx+1,nx2+1,nx3+1,nx4+1);   % at cell edge
X4 = zeros(nx+1,nx2+1,nx3+1,nx4+1);   % at cell edge
data0cc = zeros(nx,nx2,nx3,nx4); % at cell center   
for m=1:numProcs

    gridMap = vecMap(offsetsMap(m)+1:offsetsMap(m+1));
    thisMapX  = gridMap(1:end/4);
    thisMapX2 = gridMap(end/4+1:end/2);
    thisMapX3 = gridMap(end/2+1:3*end/4);
    thisMapX4 = gridMap(3*end/4+1:end);


    %%%   formulate grid at cell edges
    %
    i0 = lo_i_Map(m);
    i1 = hi_i_Map(m)+1;
    nXsub = i1-i0+1;
    j0 = lo_j_Map(m);
    j1 = hi_j_Map(m)+1;
    nX2sub = j1-j0+1;
    k0 = lo_k_Map(m);
    k1 = hi_k_Map(m)+1;
    nX3sub = k1-k0+1;
    l0 = lo_l_Map(m);
    l1 = hi_l_Map(m)+1;
    nX4sub = l1-l0+1;
    X(i0:i1,j0:j1,k0:k1,l0:l1) = reshape(thisMapX,nXsub,nX2sub,nX3sub,nX4sub);
    X2(i0:i1,j0:j1,k0:k1,l0:l1) = reshape(thisMapX2,nXsub,nX2sub,nX3sub,nX4sub);
    X3(i0:i1,j0:j1,k0:k1,l0:l1) = reshape(thisMapX3,nXsub,nX2sub,nX3sub,nX4sub);
    X4(i0:i1,j0:j1,k0:k1,l0:l1) = reshape(thisMapX4,nXsub,nX2sub,nX3sub,nX4sub);


    %%%   formulate function matrix on cell-center grid
    %
    i0data = offsets(m)+1;
    i1data = offsets(m+1);
    i0 = lo_i(m);
    i1 = hi_i(m);
    nXsub = i1-i0+1;
    j0 = lo_j(m);
    j1 = hi_j(m);
    nX2sub = j1-j0+1;
    k0 = lo_k(m);
    k1 = hi_k(m);
    nX3sub = k1-k0+1;
    l0 = lo_l(m);
    l1 = hi_l(m);
    nX4sub = l1-l0+1;
    data0cc(i0:i1,j0:j1,k0:k1,l0:l1) = reshape(vecData(i0data:i1data),nXsub,nX2sub,nX3sub,nX4sub);
    map0cc(i0:i1,j0:j1,k0:k1,l0:l1) = ones(size(data0cc(i0:i1,j0:j1,k0:k1,l0:l1)));

end


%%%   use binary map0cc to determine indices for different blocks
%
nx = length(data0cc(:,1,1,1));
ny = length(data0cc(1,:,1,1));
nvpar = length(data0cc(1,1,:,1));
nmu   = length(data0cc(1,1,1,:));
thisBox = 1;
for i=1:nx
    for j=1:ny
        if(map0cc(i,j,1,1)==1) % get lower indices for thisBox
            i0(thisBox) = i;
            j0(thisBox) = j;
            k0(thisBox) = 1;
            l0(thisBox) = 1;
            %
            for j2=j:ny    % get upper y index for thisBox
                if(map0cc(i,j2,1,1)==0)
                    j1(thisBox) = j2;
                    break;
                end
                if(j2==ny)  
                    j1(thisBox) = j2+1;
                    break;
                end
            end
            for i2=i:nx    % get upper x index for thisBox
                if(map0cc(i2,j,1,1)==0)
                    i1(thisBox) = i2;
                    break;
                end
                if(i2==nx)
                    i1(thisBox) = i2+1;
                    break;
                end
            end    
            k1(thisBox) = nvpar+1;
            l1(thisBox) = nmu+1;

            %%%   zero out thisBox so don't find it twice
            %
            i00 = i0(thisBox); i11 = i1(thisBox);
            j00 = j0(thisBox); j11 = j1(thisBox);
            k00 = k0(thisBox); k11 = k1(thisBox);
            l00 = l0(thisBox); l11 = l1(thisBox);
            map0cc(i00:i11,j00:j11,k00:k11,l00:l11) = zeros(i11-i00+1,j11-j00+1,k11-k00+1,l11-l00+1);
            thisBox = thisBox+1;
        end
    end
end
 
%%%  write grid and data in output structure for each block
%
numBlocks  = length(i0); % number of blocks
for b=1:numBlocks
    nxb = i1(b)-i0(b);
    nx2b = j1(b)-j0(b);
    nx3b = k1(b)-k0(b);
    nx4b = l1(b)-l0(b);
    %
    %
    DataStr.block(b).Xce = X(i0(b):i1(b),j0(b):j1(b),:,:);
    DataStr.block(b).X2ce = X2(i0(b):i1(b),j0(b):j1(b),:,:);
    DataStr.block(b).X3ce = X3(i0(b):i1(b),j0(b):j1(b),:,:);
    DataStr.block(b).X4ce = X4(i0(b):i1(b),j0(b):j1(b),:,:);
    DataStr.block(b).i0 = i0(b);
    DataStr.block(b).i1 = i1(b);
    DataStr.block(b).j0 = j0(b);
    DataStr.block(b).j1 = j1(b);
    DataStr.block(b).k0 = k0(b);
    DataStr.block(b).k1 = k1(b);
    DataStr.block(b).l0 = l0(b);
    DataStr.block(b).l1 = l1(b);
 %   DataStr.block(b).Fce = zeros(nxb+1,nx2b+1,nx3b+1,nx4b+1);

    %%% extend cc data by one in each dim
    %%% for compatability with ce grid
    %
    data0ce = zeros(nxb+1,nx2b+1,nx3b+1,nx4b+1);
    data0ce(1:end-1,1:end-1,1:end-1,1:end-1) = data0cc(i0(b):i1(b)-1,j0(b):j1(b)-1,k0(b):k1(b)-1,l0(b):l1(b)-1);
    data0ce(end,:,:,:) = data0ce(end-1,:,:,:);
    data0ce(:,end,:,:) = data0ce(:,end-1,:,:);
    data0ce(:,:,end,:) = data0ce(:,:,end-1,:);
    DataStr.block(b).Fce = data0ce;

end


%%%   write output information to a data structure
%
DataStr.time = time;     % time vector
DataStr.Xce = X;         % X-grid at cell edges
DataStr.X2ce = X2;         % X2-grid at cell edges
DataStr.X3ce = X3;         % X3-grid at cell edges
DataStr.X4ce = X4;         % X4-grid at cell edges
DataStr.Fcc = data0cc;     % data values at cell-center
DataStr.numProcs = numProcs;   % number of processors
DataStr.numBlocks = numBlocks;
DataStr.numGhostX = ghost.intvecti;
DataStr.numGhostY = ghost.intvectj;

end


