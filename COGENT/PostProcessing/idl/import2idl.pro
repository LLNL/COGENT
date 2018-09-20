pro import2idl, p, path=path, filename=filename
;
;
;

if not KEYWORD_SET(path) then path="./"
if not KEYWORD_SET(filename) then filename="plt.1.hydrogen.density0002.2d.hdf5"




;;-data file
fname_d=filename

;;-map file
strarr=STRSPLIT(filename, /ext, '.')
nstr=N_ELEMENTS(strarr)
fname_m=STRJOIN([strarr[0:nstr-2],'map',strarr[nstr-1]],'.')


print, "Using path ", path
print, "Importing data from ", fname_d
print, "Importing data from ", fname_m


;;-data
filename=path+fname_d
fid = H5f_OPEN(filename)
resd = H5_PARSE(filename,/read)
H5F_CLOSE, fid

;;-map
filename=path+fname_m
fid = H5f_OPEN(filename)
resm = H5_PARSE(filename,/read)
H5F_CLOSE, fid


prob_domain = resd.level_0.prob_domain._data
ghost = resd.level_0.data_attributes.ghost._data
outputGhost = resd.level_0.data_attributes.outputGhost._data
nx = prob_domain.hi_i-prob_domain.lo_i+1
ny = prob_domain.hi_j-prob_domain.lo_j+1
procs = resd.level_0.processors._data
numProcs = N_ELEMENTS(procs)

vecData = resd.level_0.data_datatype_0._data
offsets = resd.level_0.data_offsets_0._data
boxes   = resd.level_0.boxes._data
lo_i=boxes.lo_i-ghost.intvecti
lo_j=boxes.lo_j-ghost.intvectj
hi_i=boxes.hi_i-ghost.intvecti
hi_j=boxes.hi_j-ghost.intvectj

vecMap = resm.level_0.data_datatype_0._data
offsetsMap = resm.level_0.data_offsets_0._data
boxesMap   = resm.level_0.boxes._data
lo_i_map=boxes.lo_i-ghost.intvecti
lo_j_map=boxes.lo_j-ghost.intvectj
hi_i_map=boxes.hi_i-ghost.intvecti
hi_j_map=boxes.hi_j-ghost.intvectj

;;-save all data in a structure of boxes
p = CREATE_STRUCT('title', fname_d)
iglob=0
iglobm=0


for ip=0,numProcs-1 do begin

   ;;-box dimension without guard cells
   nxsub=hi_i_map[ip]-lo_i_map[ip]+1
   nysub=hi_j_map[ip]-lo_j_map[ip]+1

   ;;-box dimension with guard cells
   nxsubg=nxsub+2*outputghost.intvecti
   nysubg=nysub+2*outputghost.intvectj

   ;;-extract data from this box
   fsub=vecData[iglob:iglob+(nxsubg)*(nysubg)-1]
   iglob=iglob+(nxsubg)*(nysubg)

   ;;-extract coords from this box
   xysub=vecMap[iglobm:iglobm+2*(nxsubg+1)*(nysubg+1)-1]
   xsub=xysub[0:(nxsubg+1)*(nysubg+1)-1]
   ysub=xysub[(nxsubg+1)*(nysubg+1):*]
   iglobm=iglobm+2*(nxsubg+1)*(nysubg+1)

   ;;-cast to 2D grid and data for this box
   x2d=REFORM(xsub,nxsubg+1, nysubg+1)
   y2d=REFORM(ysub,nxsubg+1, nysubg+1)
   f2d=REFORM(fsub,nxsubg, nysubg)

   ;;-if guard cells are included they should be removed here
   ixmin=outputghost.intvecti
   ixmax=nxsubg-outputghost.intvecti
   iymin=outputghost.intvectj
   iymax=nysubg-outputghost.intvectj

   x2d=x2d[ixmin:ixmax,iymin:iymax]
   y2d=y2d[ixmin:ixmax,iymin:iymax]
   f2d=f2d[ixmin:ixmax-1,iymin:iymax-1]

   box={x2d:x2d, y2d:y2d, f2d:f2d}

   ;;-add this box to the data structure of boxes
   p = CREATE_STRUCT(p, 'box'+STRTRIM(STRING(ip),2), box)

   ;;;STOP
endfor

;
;
;
if keyword_set(DEBUG) then STOP
end
