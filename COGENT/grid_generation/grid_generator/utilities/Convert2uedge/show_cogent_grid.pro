pro show_points, domain, wexternal=wexternal, color=color
;
;
;

npsi_ext = domain.npsi_ext
npsi     = domain.npsi
npsi_tot = npsi+2*npsi_ext

nthe_ext = domain.nthe_ext
nthe     = domain.nthe
nthe_tot = nthe+2*nthe_ext



if keyword_set(WEXTERNAL) then begin
    imin=0
    imax=npsi_tot-1
    jmin=0
    jmax=nthe_tot-1
endif else begin
    imin=npsi_ext
    imax=npsi+npsi_ext-1
    jmin=nthe_ext
    jmax=nthe+nthe_ext-1
endelse


r=domain.r2[imin:imax, jmin:jmax]
z=domain.z2[imin:imax, jmin:jmax]


;;-show grid points
oplot, r, z, psym=4, syms=0.5, color=color


;
;
;
end





pro show_cells, cells, wexternal=wexternal, wtime=wtime, color=color, fill=fill, debug=debug
;
;
;

;;-number of guard cells
npsi_ext=cells.npsi_ext              
nthe_ext=cells.nthe_ext

;;-number of all cells
npsi_tot=cells.npsi+2*cells.npsi_ext 
nthe_tot=cells.nthe+2*cells.nthe_ext




if keyword_set(WEXTERNAL) then begin
    imin=0
    imax=npsi_tot-1
    jmin=0
    jmax=nthe_tot-1
endif else begin
    imin=npsi_ext
    imax=npsi_tot-npsi_ext-1
    jmin=nthe_ext
    jmax=nthe_tot-nthe_ext-1
endelse


for i0=imin,imax do begin
    for j0=jmin,jmax do begin

        if keyword_set(FILL) then begin
            POLYFILL, cells.r[i0,j0,[1,2,4,3,1]], cells.z[i0,j0,[1,2,4,3,1]], color=color
        endif else begin
            oplot, cells.r[i0,j0,[1,2,4,3,1]], cells.z[i0,j0,[1,2,4,3,1]], color=color
        endelse

        if keyword_set(WTIME) then WAIT, wtime
        ;;STOP
    endfor
    ;;print, "i0=", i0
    if keyword_set(WTIME) then WAIT, wtime
endfor


;
;
;
end



pro show_cogent_grid, alldomain, stop=stop, points=points, fill=fill, wexternal=wexternal, wtime=wtime, xrange=xrange, yrange=yrange, debug=debug
;
;
;

;;if not keyword_set(WTIME) then wtime=1e-10

ndomain = N_ELEMENTS(alldomain)
print, "Total number of domains is: ", ndomain
tek_color



;;!P.MULTI=[0,4,2,0,0]

if not keyword_set(XRANGE) then xrange=[0.6,2.5]
if not keyword_set(YRANGE) then yrange=[-1.8,1.1]

plot, xrange, yrange, /nod,/iso,/xst,/yst



for i=0,ndomain-1 do begin

    ;;-extract one domain from the list
    d = *alldomain[i]
    ;STOP
    ;;-calculate grid cells for this domain
    Domain_Cells, d, c

    ;;-show cells for this domain
    Show_Cells, c, wexternal=wexternal, wtime=wtime, fill=fill, col=i+2

    ;;-show grid points
    if keyword_set(POINTS) then Show_points, d
   

    if keyword_set(STOP) then STOP

endfor

if keyword_set(DEBUG) then STOP
;
;
;
end
