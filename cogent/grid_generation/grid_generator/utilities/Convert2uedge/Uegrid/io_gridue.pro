pro read_gridue, g
;
;
;

;;print, "In read_gridue (IDL)..."

nxm=0l
nym=0l
status = CALL_EXTERNAL('io_gridue.so', 'read_gridue_dims', nxm,nym)
print, "nxm,nym=", nxm, nym
;;;stop


ixpt1=0l
ixpt2=0l
iysptrx1=0l

rm=dblarr(nxm+2,nym+2,5)
zm=dblarr(nxm+2,nym+2,5)
psi=dblarr(nxm+2,nym+2,5)
br=dblarr(nxm+2,nym+2,5)
bz=dblarr(nxm+2,nym+2,5)
bpol=dblarr(nxm+2,nym+2,5)
bphi=dblarr(nxm+2,nym+2,5)
b=dblarr(nxm+2,nym+2,5)


status = CALL_EXTERNAL('io_gridue.so', 'read_gridue', $
                        nxm,nym,ixpt1,ixpt2,iysptrx1, rm,zm,psi,br,bz,bpol,bphi,b)
;;stop
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
    bpol:bpol,$
    bphi:bphi,$
    b:b}

;
;
;
end




pro write_gridue, g
;
;
;

nxm=LONG(g.nxm)
nym=LONG(g.nym)
ixpt1=LONG(g.ixpt1)
ixpt2=LONG(g.ixpt2)
iysptrx1=LONG(g.iysptrx1)
rm=DOUBLE(g.rm)
zm=DOUBLE(g.zm)
psi=DOUBLE(g.psi)
br=DOUBLE(g.br)
bz=DOUBLE(g.bz)
bpol=DOUBLE(g.bpol)
bphi=DOUBLE(g.bphi)
b=DOUBLE(g.b)


;;print, "In write_gridue.pro: rm[0,0,0]=", rm[0,0,0]

status = CALL_EXTERNAL('io_gridue.so', 'write_gridue', $
                       nxm,nym,ixpt1,ixpt2,iysptrx1, rm,zm,psi,br,bz,bpol,bphi,b)


;
;
;
end


pro plot_gridue, g, xrange=xrange, yrange=yrange, debug=debug
;
;
;

plot, g.rm, g.zm, xrange=xrange, yrange=yrange, /iso, /nodata

for i=0,g.nxm+1 do begin
   for j=0,g.nym+1 do begin

      oplot, g.rm[i,j,[1,2,4,3,1]], g.zm[i,j,[1,2,4,3,1]]

   endfor
endfor

for i=0,0 do begin
   for j=0,g.nym+1 do begin
      oplot, g.rm[i,j,[1,2,4,3,1]], g.zm[i,j,[1,2,4,3,1]], col=2
   endfor
endfor

for i=g.nxm+1,g.nxm+1 do begin
   for j=0,g.nym+1 do begin
      oplot, g.rm[i,j,[1,2,4,3,1]], g.zm[i,j,[1,2,4,3,1]], col=2
   endfor
endfor


;;-show radial guard cells
for i=0,g.nxm+1 do begin
   for j=0,0 do begin

      oplot, g.rm[i,j,[1,2,4,3,1]], g.zm[i,j,[1,2,4,3,1]], col=3

   endfor
endfor

for i=0,g.nxm+1 do begin
   for j=g.nym+1,g.nym+1 do begin

      oplot, g.rm[i,j,[1,2,4,3,1]], g.zm[i,j,[1,2,4,3,1]], col=3

   endfor
endfor

;
;
;
if keyword_set(DEBUG) then STOP
end
