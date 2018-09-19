pro plot_subgrid, b, color=color
;

  if NOT(KEYWORD_SET(color)) then color=1

  sz=SIZE(b.f2d)
  nx=sz[1]
  ny=sz[2]

  for i=0,nx-1 do begin
     for j=0,ny-1 do begin
        x=[b.x2d[i,j],b.x2d[i+1,j],b.x2d[i+1,j+1],b.x2d[i,j+1],b.x2d[i,j]]
        y=[b.y2d[i,j],b.y2d[i+1,j],b.y2d[i+1,j+1],b.y2d[i,j+1],b.y2d[i,j]]
        oplot, x,y, color=color
     endfor
  endfor

;
;
end


pro show_grid, p, xrange=xrange, yrange=yrange, wcolor=wcolor, debug=debug
;
;
;

TEK_COLOR


tags=TAG_NAMES(p)
ntags=N_ELEMENTS(tags)



;;-determine the range of grid
for i=1,ntags-1 do begin
   box=p.(i)
   if (i eq 1) then begin
      xmin=MIN(box.x2d)
      xmax=MAX(box.x2d)
      ymin=MIN(box.y2d)
      ymax=MAX(box.y2d)
   endif else begin
      xmin=MIN([xmin,MIN(box.x2d)])
      xmax=MAX([xmax,MAX(box.x2d)])
      ymin=MIN([ymin,MIN(box.y2d)])
      ymax=MAX([ymax,MAX(box.y2d)])
   endelse
endfor


if not keyword_set(xrange) then xrange=[xmin,xmax]
if not keyword_set(yrange) then yrange=[ymin,ymax]
plot, xrange, yrange, /nod,/xst,/yst,/iso
XYOUTS,/norm,0,0, p.title




for i=1,ntags-1 do begin

   box=p.(i)

   if KEYWORD_SET(wcolor) then begin
      Plot_subgrid, box, color=2+(i mod 8)
   endif else begin
      Plot_subgrid, box
   endelse

endfor


;
;
;
if keyword_set(DEBUG) then STOP
end
