pro remap_domain, indomain=indomain, outdomain=outdomain, npsi_new=npsi_new, nthe_new=nthe_new, debug=debug
;
; - Remamp a domain grid points to new indices but within old physical
;   boundaries
;
;======================================================================;


npsi = indomain.npsi
npsi_ext = indomain.npsi_ext
nthe = indomain.nthe
nthe_ext = indomain.nthe_ext


;;-new number of internal grid points
if not keyword_set(NPSI_NEW) then npsi_new=npsi
if not keyword_set(NTHE_NEW) then nthe_new=nthe


;;-2D array of all original grid points (including guard points)
r2 = indomain.r2
z2 = indomain.z2
rBr2 = indomain.rBr2
rBz2 = indomain.rBz2
psi2 = indomain.psi2


;;-interpolate to a new set of grid points (including guard points)
ipsi_new = [findgen(npsi_ext), npsi_ext+(npsi-1)*FINDGEN(npsi_new)/(npsi_new-1), npsi_ext+npsi+findgen(npsi_ext)]
ithe_new = [findgen(nthe_ext), nthe_ext+(nthe-1)*FINDGEN(nthe_new)/(nthe_new-1), nthe_ext+nthe+findgen(nthe_ext)]


r2new = INTERPOLATE(r2, ipsi_new, ithe_new, /grid)
z2new = INTERPOLATE(z2, ipsi_new, ithe_new, /grid)
rBr2new = INTERPOLATE(rBr2, ipsi_new, ithe_new, /grid)
rBz2new = INTERPOLATE(rBz2, ipsi_new, ithe_new, /grid)
psi2new = INTERPOLATE(psi2, ipsi_new, ithe_new, /grid)


outdomain={npsi_ext:npsi_ext, nthe_ext:nthe_ext, npsi:npsi_new, nthe:nthe_new, $
           r2:r2new, z2:z2new, rBr2:rBr2new, rBz2:rBz2new, psi2:psi2new}


;
;
;
if keyword_set(debug) then STOP
end




pro remap_cogent_grid, alldomain, alldomain2
;
; -Remap the PF region to make the number of radial grid cells match
; that in the core, like in UEDGE
;===================================================================;;

alldomain2=alldomain

d0 = *alldomain[0]
ncore = d0.npsi

d6 = *alldomain[6]
Remap_domain, in=d6, out=d6new, npsi_new=ncore
alldomain2[6] = PTR_NEW(d6new)

d7 = *alldomain[7]
Remap_domain, in=d7, out=d7new, npsi_new=ncore
alldomain2[7] = PTR_NEW(d7new)




;
;
;
end
