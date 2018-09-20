function [DataStr] = import2Ddata(dataPath,species)


if nargin<2 
  species = 1; % if species not specified, use species=1
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   import 2D data from COGENT output folder
%%%
%%%   input parameters:
%%%   dataPath = 'string/to/plt_something_plots/' 
%%%   species = species number (default if not specified is 1 
%%%                             and use 0 for non-species outputs)
%%%
%%%   output parameter is a structure with:
%%%   block structure with Xce, Yce, and Fce(time) for each block
%%%   time array
%%%   cell-edge grid (Xce and Yce) for all blocks
%%%   cell-center data (Fcc(time)) for all blocks
%%%   number of components
%%%   number of processors
%%%   number of blocks
%%%   number of ghost cells at each end in X
%%%   number of ghost cells at each end in Y
%%%
%%%   block structure has local Xce, Yce, and Fce(time) for each block, where 
%%%   Fce is cell centered data extened by one in each direction
%%%   to be compatible with Xce and Yce for pcolor and contour plots
%%%
%%%   Note that the function value is defined at cell-center
%%%   while the grid is defined at cell-edge
%%%
%%%   April 25, 2018
%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(species==0) % non species outputs
    plisting    = dir([dataPath,'plt.*.2d.hdf5']);
    plistingMap = dir([dataPath,'plt.*.2d.map.hdf5']);
else
    plisting    = dir([dataPath,'*.',num2str(species),'*.2d.hdf5']);
    plistingMap = dir([dataPath,'*.',num2str(species),'*.2d.map.hdf5']);
end

if(isempty(plisting))
    error(['Error: could not find hdf5 files for ' ...
           'species = ', num2str(species),' in dataPath = ',dataPath]);
end


%%%   reorder filenames to be time sequential (not always the case
%%%   when using plisting)
%
time = zeros(1,length(plisting));
for i=1:length(plisting)
    thisFile = [dataPath,plisting(i).name];
    time(i) = h5readatt(thisFile,'/level_0','time');
    if(i==length(plisting) ) % get grid data
      %  dx = h5readatt(thisFile,'/level_0','dx');
        prob_domain = h5readatt(thisFile,'/level_0','prob_domain');
        ghost = h5readatt(thisFile,'/level_0/data_attributes','ghost');
        outputGhost = h5readatt(thisFile,'/level_0/data_attributes','outputGhost');
        %nDims = length(fileinfo.GroupHierarchy.Groups(2).Attributes(5).Value.Data)/2;
        nComps = h5readatt(thisFile,'/','num_components');
        nx = double(prob_domain.hi_i-prob_domain.lo_i+1);
        ny = double(prob_domain.hi_j-prob_domain.lo_j+1);
        procs = hdf5read(thisFile,'/level_0/Processors');
        numProcs = length(procs);
    end
end
[time,index] = sort(time);
plisting = plisting(index);
plistingMap = plistingMap(index);


%%%   loop over time and load the data
%
nt = length(plisting);
data0out = zeros(nx+2*ghost.intvecti,ny+2*ghost.intvectj,nComps,nt);
for it=1:nt

    thisFile = [dataPath,plisting(it).name];
    fileinfo = hdf5info(thisFile);
    vecData  = hdf5read(thisFile,'/level_0/data:datatype=0');
    offsets  = hdf5read(thisFile,'/level_0/data:offsets=0');
    boxes    = hdf5read(thisFile,'/level_0/boxes');
    for iP = 1:numProcs
        lo_i(iP) = boxes(iP).Data{1}+1-ghost.intvecti;
        lo_j(iP) = boxes(iP).Data{2}+1-ghost.intvectj;       
        hi_i(iP) = boxes(iP).Data{3}+1+ghost.intvecti;
        hi_j(iP) = boxes(iP).Data{4}+1+ghost.intvectj;
    end
    
    %%%   some indices are negative, so shift for MATLAB
    %
    min_lo_i = min(lo_i);
    lo_i = lo_i-double(min_lo_i)+1;
    hi_i = hi_i-double(min_lo_i)+1;
    min_lo_j = min(lo_j);
    lo_j = lo_j-double(min_lo_j)+1;
    hi_j = hi_j-double(min_lo_j)+1;
    
    
    thisFileMap = [dataPath,plistingMap(it).name];
    fileinfoMap = hdf5info(thisFileMap);
    vecMap      = hdf5read(thisFileMap,'/level_0/data:datatype=0');
    offsetsMap  = hdf5read(thisFileMap,'/level_0/data:offsets=0');
    boxesMap    = hdf5read(thisFileMap,'/level_0/boxes');
    lo_i_Map = zeros(1,numProcs);
    lo_j_Map = zeros(1,numProcs);
    hi_i_Map = zeros(1,numProcs);
    hi_j_Map = zeros(1,numProcs);
    for iP = 1:numProcs
        lo_i_Map(iP) = boxesMap(iP).Data{1}+1-ghost.intvecti;
        lo_j_Map(iP) = boxesMap(iP).Data{2}+1-ghost.intvectj;       
        hi_i_Map(iP) = boxesMap(iP).Data{3}+1+ghost.intvecti;
        hi_j_Map(iP) = boxesMap(iP).Data{4}+1+ghost.intvectj;
    end
    min_lo_i_Map = min(lo_i_Map);
    lo_i_Map = lo_i_Map-double(min_lo_i_Map)+1;
    hi_i_Map = hi_i_Map-double(min_lo_i_Map)+1;
    min_lo_j_Map = min(lo_j_Map);
    lo_j_Map = lo_j_Map-double(min_lo_j_Map)+1;
    hi_j_Map = hi_j_Map-double(min_lo_j_Map)+1;
    
    
    %%%   Note that boxes and boxesMap are always the same, while   %%%
    %%%   data is at cell center and grid at cell-edges             %%%
    

    %%%   map to reshaped grid
    %
    X = zeros(nx+1,ny+1);   % at cell edge
    Y = zeros(nx+1,ny+1);   % at cell edge
    data0cc = zeros(nx,ny); % at cell center   
    for m=1:numProcs

        gridMap = vecMap(offsetsMap(m)+1:offsetsMap(m+1));
        thisMapX = gridMap(1:end/2);
        thisMapY = gridMap(end/2+1:end);
        

        %%%   formulate grid at cell edges
        %
        i0 = lo_i_Map(m);
        i1 = hi_i_Map(m)+1;
        nXsub = i1-i0+1;
        j0 = lo_j_Map(m);
        j1 = hi_j_Map(m)+1;
        nYsub = j1-j0+1;
        X(i0:i1,j0:j1) = reshape(thisMapX,nXsub,nYsub);
        Y(i0:i1,j0:j1) = reshape(thisMapY,nXsub,nYsub);

        
        %%%   formulate function matrix on cell-center grid
        %
        i0data = offsets(m)+1;
        i1data = offsets(m+1);
        i0 = lo_i(m);
        i1 = hi_i(m);
        nXsub = i1-i0+1;
        j0 = lo_j(m);
        j1 = hi_j(m);
        nYsub = j1-j0+1;
        
        thisvecData = vecData(i0data:i1data);
        subSize = length(thisvecData)/nComps;
        for n=1:nComps
            nlow = 1+(n-1)*subSize;
            nup = n*subSize;
            thisvecData0 = thisvecData(nlow:nup);
            data0cc(i0:i1,j0:j1,n) = reshape(thisvecData0,nXsub,nYsub);       
        end
        map0cc(i0:i1,j0:j1) = ones(size(squeeze(data0cc(i0:i1,j0:j1,1))));   
        
        %data0cc(i0:i1,j0:j1) = reshape(vecData(i0data:i1data),nXsub,nYsub);
        %map0cc(i0:i1,j0:j1) = ones(size(data0cc(i0:i1,j0:j1)));
        
    end

    if(it==1)
    %%%   use binary map0cc to determine indices for different blocks
    %
    nx = length(data0cc(:,1,1));
    ny = length(data0cc(1,:,1));
    thisBox = 1;
    for i=1:nx
        for j=1:ny
            if(map0cc(i,j)==1) % get lower indices for thisBox
                i0(thisBox) = i;
                j0(thisBox) = j;
                %
                for j2=j:ny    % get upper y index for thisBox
                    if(map0cc(i,j2)==0)
                        j1(thisBox) = j2;
                        break;
                    end
                    if(j2==ny)  
                        j1(thisBox) = j2+1;
                        break;
                    end
                end
                for i2=i:nx    % get upper x index for thisBox
                    if(map0cc(i2,j)==0)
                        i1(thisBox) = i2;
                        break;
                    end
                    if(i2==nx)
                        i1(thisBox) = i2+1;
                        break;
                    end
                end    

                %%%   zero out thisBox so don't find it twice
                %
                i00 = i0(thisBox); i11 = i1(thisBox);
                j00 = j0(thisBox); j11 = j1(thisBox);
                map0cc(i00:i11,j00:j11) = zeros(i11-i00+1,j11-j00+1);
                thisBox = thisBox+1;
            end
        end
    end
    
    %%%  write grid and data in output structure for each block
    %
    numBlocks  = length(i0); % number of blocks
    for b=1:numBlocks
        nxb = i1(b)-i0(b);
        nyb = j1(b)-j0(b);
        %
        %
        DataStr.block(b).Xce = X(i0(b):i1(b),j0(b):j1(b));
        DataStr.block(b).Yce = Y(i0(b):i1(b),j0(b):j1(b));
        DataStr.block(b).i0 = i0(b);
        DataStr.block(b).i1 = i1(b);
        DataStr.block(b).j0 = j0(b);
        DataStr.block(b).j1 = j1(b);      
        DataStr.block(b).Fce = squeeze(zeros(nxb+1,nyb+1,nComps,nt));
    end
    
    end % end if nt==1 statement
    
    for b=1:numBlocks
        
        %%% extend cc data by one in each dim
        %%% for compatability with ce grid
        %
        i0 = DataStr.block(b).i0;
        i1 = DataStr.block(b).i1;
        j0 = DataStr.block(b).j0;
        j1 = DataStr.block(b).j1;
        nxb = i1-i0;
        nyb = j1-j0;
        %
        data0ce = zeros(nxb+1,nyb+1,nComps);
        data0ce(1:end-1,1:end-1,:) = data0cc(i0:i1-1,j0:j1-1,:);
        data0ce(end,:,:) = data0ce(end-1,:,:);
        data0ce(:,end,:) = data0ce(:,end-1,:);
        if(nComps==1)
            DataStr.block(b).Fce(:,:,it) = squeeze(data0ce);
        else
            DataStr.block(b).Fce(:,:,:,it) = data0ce;
        end
        
    end
  %  size(data0cc)
  %  size(data0out)
    data0out(:,:,:,it) = data0cc;
    
end % end time loop


%%%   write output information to a data structure
%
DataStr.time = time;     % time vector
DataStr.Xce = X;         % X-grid at cell edges
DataStr.Yce = Y;         % Y-grid at cell edges
DataStr.Fcc = squeeze(data0out);  % data values at cell-center
DataStr.nComps = nComps;       % number of components for data
DataStr.numProcs = numProcs;   % number of processors
DataStr.numBlocks = numBlocks; % number of blocks
DataStr.numGhostX = ghost.intvecti;
DataStr.numGhostY = ghost.intvectj;

end


