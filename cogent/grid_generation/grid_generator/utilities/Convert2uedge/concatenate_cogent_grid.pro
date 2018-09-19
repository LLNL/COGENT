function nan_present, arr
;
;;-return 1 if array contains a nan on infinity, 0 otherwise
;===========================================================;;

nan = WHERE(FINITE(arr) ne 1)

if ((nan[0] ne -1)) then begin
    res=1
endif else begin 
    res=0
endelse

return, res
end




pro domain_patch, domain, r=r, z=z, brr=brr, bzr=bzr, psi=psi, guardpsi=guardpsi, guardthe=guardthe, debug=debug
;
;;-extract a subset of cells, with one guard cells include on external
;;boundaries, transpose to match the UEDGE convention for poloidal
;;index (increases going clockwise around the plasma)
;;
;;-Note: choosing keyword names to make them unambigious
;;====================================================================;;;

;;print, guardthe

Domain_cells, domain, c



case guardpsi of

    'top': begin
        ipsi = c.npsi_ext + INDGEN(c.npsi+1)
    end

    'bot': begin
        ;;-include guards on the top boundary
        ipsi = c.npsi_ext-1 + INDGEN(c.npsi+1)
    end

    else: begin
        ipsi = c.npsi_ext + INDGEN(c.npsi)        
    end

endcase



case guardthe of

    'top': begin
        print, 'the top'
        ithe = REVERSE(c.nthe_ext-1 + INDGEN(c.nthe+1))
    end

    'bot': begin
        print, 'the bot'
        ithe = REVERSE(c.nthe_ext + INDGEN(c.nthe+1))
    end

    else: begin
        print, 'the none'
        ithe = REVERSE(c.nthe_ext + INDGEN(c.nthe))
    end

endcase



r=TRANSPOSE(c.r[ipsi,ithe,*],[1,0,2])
z=TRANSPOSE(c.z[ipsi,ithe,*],[1,0,2])
Brr=TRANSPOSE(c.rBr[ipsi,ithe,*],[1,0,2])
Bzr=TRANSPOSE(c.rBz[ipsi,ithe,*],[1,0,2])
psi=TRANSPOSE(c.psi[ipsi,ithe,*],[1,0,2])


if (Nan_Present(r) or Nan_Present(z) or Nan_Present(brr) or $
    Nan_Present(bzr) or Nan_Present(psi) ) then begin
    print, "Detected NaN or Infinity in data, stopping..."
    STOP
endif

;
;
end






pro concatenate_cogent_grid, alldomain, g, RBTOR=RBTOR, export=export, debug=debug
;
;
;

if not keyword_set(RBTOR) then begin
    rBtor =-3.5 ;; Tesla*m
    Result = DIALOG_MESSAGE( 'Using RxBtor='+STRING(rBtor, format='(f10.4)')+' [mxT]')
endif


;;-number of subdomains
ndomain = N_ELEMENTS(alldomain)

if (ndomain ne 8) then begin
    print, "Expecting 8 subdomains, given ", ndomain
    STOP
endif


;;-total number of internal cell in psi and theta
nym = 0 ;;-radial
nxm = 0 ;;-poloidal


for i=0,7 do begin

    domain = *alldomain[i]

    if (i eq 4 or i eq 5) then begin
        ;;- include one poloidal guard cell
        nxm =  nxm + domain.nthe
    end

    if (i eq 2 or i eq 3) then begin
        ;;- include only internal cells
        nxm =  nxm + domain.nthe-1
    end    

    if (i eq 0 or i eq 3) then begin
        ;;- include one radial guard cell
        nym =  nym + domain.npsi
    end

endfor


;;-UEDGE conventions
nxm=nxm-2
nym=nym-2


print, "The total domain dimensions"
print, "nym=", nym+2
print, "nxm=", nxm+2

;;-structure containing the whole LSN domain
rm=dblarr(nxm+2,nym+2,5)
zm=dblarr(nxm+2,nym+2,5)
psi=dblarr(nxm+2,nym+2,5)
br=dblarr(nxm+2,nym+2,5)
bz=dblarr(nxm+2,nym+2,5)
rbr=dblarr(nxm+2,nym+2,5)
rbz=dblarr(nxm+2,nym+2,5)
bpol=dblarr(nxm+2,nym+2,5)
bphi=dblarr(nxm+2,nym+2,5)
b=dblarr(nxm+2,nym+2,5)
ixpt1=1
ixpt2=1
iysptrx1=1





;;-put data from subdomains outside of the separatrix
iysptrx1=(*alldomain[6]).npsi-1
ixpt1=(*alldomain[4]).nthe-1
ixpt2=(*alldomain[2]).nthe+(*alldomain[3]).nthe+(*alldomain[4]).nthe-3


g={$
    nxm:nxm,$
    nym:nym,$
    ixpt1:ixpt1,$
    ixpt2:ixpt2,$
    iysptrx1:iysptrx1,$
    rm:rm,$
    zm:zm,$
    psi:psi,$
    br:br,$
    bz:bz,$
    rbr:rbr,$
    rbz:rbz,$
    bpol:bpol,$
    bphi:bphi,$
    b:b}



;;STOP
;;-extract a subdomain and calculate grid cells
d=*alldomain[6]
Domain_patch, d, r=rr, z=zz, brr=rbr, bzr=rbz, psi=psi, guardpsi='bot', guardthe='bot'
g.rm[0:ixpt1,0:iysptrx1,*] = rr
g.zm[0:ixpt1,0:iysptrx1,*] = zz
g.rbr[0:ixpt1,0:iysptrx1,*] = rbr
g.rbz[0:ixpt1,0:iysptrx1,*] = rbz
g.psi[0:ixpt1,0:iysptrx1,*] = psi



d=*alldomain[4]
Domain_patch, d, r=rr, z=zz, brr=rbr, bzr=rbz, psi=psi, guardpsi='top', guardthe='bot'
g.rm[0:ixpt1,iysptrx1+1:*,*] = rr
g.zm[0:ixpt1,iysptrx1+1:*,*] = zz
g.rbr[0:ixpt1,iysptrx1+1:*,*] = rbr
g.rbz[0:ixpt1,iysptrx1+1:*,*] = rbz
g.psi[0:ixpt1,iysptrx1+1:*,*] = psi




d=*alldomain[0]
Domain_patch, d, r=rr1, z=zz1, brr=rbr1, bzr=rbz1, psi=psi1, guardpsi='bot', guardthe='none'
d=*alldomain[1]
Domain_patch, d, r=rr2, z=zz2, brr=rbr2, bzr=rbz2, psi=psi2, guardpsi='bot', guardthe='none'
g.rm[ixpt1+1:ixpt2,0:iysptrx1,*] = [rr1,rr2]
g.zm[ixpt1+1:ixpt2,0:iysptrx1,*] = [zz1,zz2]
g.rbr[ixpt1+1:ixpt2,0:iysptrx1,*] = [rbr1,rbr2]
g.rbz[ixpt1+1:ixpt2,0:iysptrx1,*] = [rbz1,rbz2]
g.psi[ixpt1+1:ixpt2,0:iysptrx1,*] = [psi1,psi2]



d=*alldomain[2]
Domain_patch, d, r=rr1, z=zz1, brr=rbr1, bzr=rbz1, psi=psi1, guardpsi='top', guardthe='none'
d=*alldomain[3]
Domain_patch, d, r=rr2, z=zz2, brr=rbr2, bzr=rbz2, psi=psi2, guardpsi='top', guardthe='none'
g.rm[ixpt1+1:ixpt2,iysptrx1+1:*,*] = [rr1,rr2]
g.zm[ixpt1+1:ixpt2,iysptrx1+1:*,*] = [zz1,zz2]
g.rbr[ixpt1+1:ixpt2,iysptrx1+1:*,*] = [rbr1,rbr2]
g.rbz[ixpt1+1:ixpt2,iysptrx1+1:*,*] = [rbz1,rbz2]
g.psi[ixpt1+1:ixpt2,iysptrx1+1:*,*] = [psi1,psi2]



d=*alldomain[7]
Domain_patch, d, r=rr, z=zz, brr=rbr, bzr=rbz, psi=psi, guardpsi='bot', guardthe='top'
g.rm[ixpt2+1:*,0:iysptrx1,*] = rr
g.zm[ixpt2+1:*,0:iysptrx1,*] = zz
g.rbr[ixpt2+1:*,0:iysptrx1,*] = rbr
g.rbz[ixpt2+1:*,0:iysptrx1,*] = rbz
g.psi[ixpt2+1:*,0:iysptrx1,*] = psi



d=*alldomain[5]
Domain_patch, d, r=rr, z=zz, brr=rbr, bzr=rbz, psi=psi, guardpsi='top', guardthe='top'
g.rm[ixpt2+1:*,iysptrx1+1:*,*] = rr
g.zm[ixpt2+1:*,iysptrx1+1:*,*] = zz
g.rbr[ixpt2+1:*,iysptrx1+1:*,*] = rbr
g.rbz[ixpt2+1:*,iysptrx1+1:*,*] = rbz
g.psi[ixpt2+1:*,iysptrx1+1:*,*] = psi




Show_Cogent_Grid, alldomain, xr=[0.8, 2.5], yr=[-2,-1]
oplot, rr[*,*,0], zz[*,*,0], psym=4
;;STOP

;;-reconstruct the missing components of the grid structure 
g.br = g.rbr/g.rm
g.bz = g.rbz/g.rm
g.bpol=SQRT(g.br^2+g.bz^2)
g.bphi = rBtor/g.rm
g.b = SQRT(g.bpol^2+g.bphi^2)


if keyword_set(EXPORT) then begin
    Write_gridue, g
endif


;
;
;
if keyword_set(DEBUG) then STOP
end
