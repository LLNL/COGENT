pro import_cogent_grid, fname=fname, alldomain
;
;
;


if not keyword_set(fname) then fname="COGENT_mapping_UE"
print, "Opening file ", fname

SPAWN, "grep nan "+fname, res
if (STRLEN(res[0]) gt 0) then begin
    print, "Found nan's in ", fname, ", stopped..."
    STOP
endif


OPENR, inunit, fname, /GET_LUN

line = ''


WHILE ~ EOF(inunit) DO BEGIN

   READF, inunit, line

   ;;-convert line to array of strings
   str= STRSPLIT(line,/ex)

   ;;-check if line contains anything other than exponential form numbers:
   res= STREGEX(str,"^([+-][0-9]+[.]?[0-9]*|[0-9]*[.][0-9]+)[eE][+-]?[0-9]+$",/boolean)

   if MAX(ABS(res-1)) then begin

       if N_ELEMENTS(domain_name) eq 0 then begin
           domain_name=line
       endif else begin
           domain_name=[domain_name,line]
       endelse

       header_line=line
       print, "Current domain header line: ", header_line

       header = STRSPLIT(str,' ',/ex)
       dname    = (header[0])[0]
       npsi     = (FIX(header[1]))[0]
       npsi_ext = (FIX(header[2]))[0]
       nthe     = (FIX(header[3]))[0]
       nthe_ext = (FIX(header[4]))[0]
       print, dname, npsi, npsi_ext, nthe, nthe_ext
       ;;STOP


       ;;-number of grid points in the current domain
       npt=(npsi+2*npsi_ext)*(nthe+2*nthe_ext)
       domain = {name:dname, $
                 r:DBLARR(npt),z:DBLARR(npt),rBr:DBLARR(npt),rBz:DBLARR(npt),psi:DBLARR(npt),$
                 r2:DBLARR(npsi+2*npsi_ext,nthe+2*nthe_ext),$
                 z2:DBLARR(npsi+2*npsi_ext,nthe+2*nthe_ext),$
                 rBr2:DBLARR(npsi+2*npsi_ext,nthe+2*nthe_ext),$
                 rBz2:DBLARR(npsi+2*npsi_ext,nthe+2*nthe_ext),$
                 psi2:DBLARR(npsi+2*npsi_ext,nthe+2*nthe_ext),$
                 npt:npt, npsi:npsi, npsi_ext:npsi_ext, nthe:nthe, nthe_ext:nthe_ext}
       
       for j=0,npt-1 do begin
           READF, inunit, r0, z0, rBr0, rBz0, psi0
           domain.r[j]=r0
           domain.z[j]=z0
           domain.rBr[j]=rBr0
           domain.rBz[j]=rBz0
           domain.psi[j]=psi0
       endfor


       ;;-cast in 2D form
       npsi_tot = npsi+2*npsi_ext
       nthe_tot = nthe+2*nthe_ext
       domain.r2   = REFORM(domain.r,   npsi_tot, nthe_tot)
       domain.z2   = REFORM(domain.z,   npsi_tot, nthe_tot)
       domain.rBr2 = REFORM(domain.rBr, npsi_tot, nthe_tot)
       domain.rBz2 = REFORM(domain.rBz, npsi_tot, nthe_tot)
       domain.psi2 = REFORM(domain.psi, npsi_tot, nthe_tot)

       ;;STOP

       ;;-add to an array of pointers to structures
       if (N_ELEMENTS(alldomain) eq 0) then begin
           alldomain = [PTR_NEW(domain)]
       endif else begin
           alldomain = [alldomain, PTR_NEW(domain)]
       endelse


   endif

ENDWHILE

 

; Close the files and deallocate the units:
FREE_LUN, inunit
CLOSE, inunit
print, "Closing file ", fname


ndomain=N_ELEMENTS(domain_name)
print, "Number of domains: ", ndomain
for i=0,ndomain-1 do begin
    print, " ", domain_name[i]
endfor


;
;
;
end




pro domain_cells, domain, cells
;
;-Calculate FV cells using UEDGE convention for a COGENT domain
;

npsi_ext = domain.npsi_ext ;;-number of external grid points
nthe_ext = domain.nthe_ext ;;-number of external grid points
npsi     = domain.npsi     ;;-number of internal grid points
nthe     = domain.nthe     ;;-number of internal grid points


cells = {npsi_ext:npsi_ext, nthe_ext:nthe_ext,$ ;;-number of guard cells
         npsi:npsi-1, nthe:nthe-1, $ ;;-number of internal cells
         r:DBLARR(npsi+2*npsi_ext-1,nthe+2*nthe_ext-1,5),$
         z:DBLARR(npsi+2*npsi_ext-1,nthe+2*nthe_ext-1,5),$
         rBr:DBLARR(npsi+2*npsi_ext-1,nthe+2*nthe_ext-1,5),$
         rBz:DBLARR(npsi+2*npsi_ext-1,nthe+2*nthe_ext-1,5),$
         psi:DBLARR(npsi+2*npsi_ext-1,nthe+2*nthe_ext-1,5)}



;;-convert from grid points to cells (using UEDGE convention)

;;-SW corner
cells.r[*,*,1] = domain.r2[0:npsi+2*npsi_ext-2, 0:nthe+2*nthe_ext-2]
cells.z[*,*,1] = domain.z2[0:npsi+2*npsi_ext-2, 0:nthe+2*nthe_ext-2]
cells.rBr[*,*,1] = domain.rBr2[0:npsi+2*npsi_ext-2, 0:nthe+2*nthe_ext-2]
cells.rBz[*,*,1] = domain.rBz2[0:npsi+2*npsi_ext-2, 0:nthe+2*nthe_ext-2]
cells.psi[*,*,1] = domain.psi2[0:npsi+2*npsi_ext-2, 0:nthe+2*nthe_ext-2]


;;-SE corner
cells.r[*,*,2] = domain.r2[0:npsi+2*npsi_ext-2, 1:nthe+2*nthe_ext-1]
cells.z[*,*,2] = domain.z2[0:npsi+2*npsi_ext-2, 1:nthe+2*nthe_ext-1]
cells.rBr[*,*,2] = domain.rBr2[0:npsi+2*npsi_ext-2, 1:nthe+2*nthe_ext-1]
cells.rBz[*,*,2] = domain.rBz2[0:npsi+2*npsi_ext-2, 1:nthe+2*nthe_ext-1]
cells.psi[*,*,2] = domain.psi2[0:npsi+2*npsi_ext-2, 1:nthe+2*nthe_ext-1]


;;-NW corner
cells.r[*,*,3] = domain.r2[1:npsi+2*npsi_ext-1, 0:nthe+2*nthe_ext-2]
cells.z[*,*,3] = domain.z2[1:npsi+2*npsi_ext-1, 0:nthe+2*nthe_ext-2]
cells.rBr[*,*,3] = domain.rBr2[1:npsi+2*npsi_ext-1, 0:nthe+2*nthe_ext-2]
cells.rBz[*,*,3] = domain.rBz2[1:npsi+2*npsi_ext-1, 0:nthe+2*nthe_ext-2]
cells.psi[*,*,3] = domain.psi2[1:npsi+2*npsi_ext-1, 0:nthe+2*nthe_ext-2]


;;-NE corner
cells.r[*,*,4] = domain.r2[1:npsi+2*npsi_ext-1, 1:nthe+2*nthe_ext-1]
cells.z[*,*,4] = domain.z2[1:npsi+2*npsi_ext-1, 1:nthe+2*nthe_ext-1]
cells.rBr[*,*,4] = domain.rBr2[1:npsi+2*npsi_ext-1, 1:nthe+2*nthe_ext-1]
cells.rBz[*,*,4] = domain.rBz2[1:npsi+2*npsi_ext-1, 1:nthe+2*nthe_ext-1]
cells.psi[*,*,4] = domain.psi2[1:npsi+2*npsi_ext-1, 1:nthe+2*nthe_ext-1]


;;-center
cells.r[*,*,0] = 0.25*TOTAL(cells.r[*,*,1:4], 3)
cells.z[*,*,0] = 0.25*TOTAL(cells.z[*,*,1:4], 3)
cells.rBr[*,*,0] = 0.25*TOTAL(cells.rBr[*,*,1:4], 3)
cells.rBz[*,*,0] = 0.25*TOTAL(cells.rBz[*,*,1:4], 3)
cells.psi[*,*,0] = 0.25*TOTAL(cells.psi[*,*,1:4], 3)

;
;
;
end
