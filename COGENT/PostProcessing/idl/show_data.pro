pro polyfill_subgrid, b, minval=minval, maxval=maxval
;

  sz=SIZE(b.f2d)
  nx=sz[1]
  ny=sz[2]

  for i=0,nx-1 do begin
     for j=0,ny-1 do begin
        x=[b.x2d[i,j],b.x2d[i+1,j],b.x2d[i+1,j+1],b.x2d[i,j+1],b.x2d[i,j]]
        y=[b.y2d[i,j],b.y2d[i+1,j],b.y2d[i+1,j+1],b.y2d[i,j+1],b.y2d[i,j]]
        ;;oplot, x,y, color=color
        color=5+249*(b.f2d[i,j]-minval)/(maxval-minval)
        ;;print, color
        ;;STOP
        POLYFILL, x,y, color=color
     endfor
  endfor

;
;
end


pro show_data, p, xrange=xrange, yrange=yrange, minval=minval, maxval=maxval, debug=debug
;
;
;

LOADCT, 39

tags=TAG_NAMES(p)
ntags=N_ELEMENTS(tags)
;;print, "Tag names: ", tags

;;plot, [0,10],[-10,10],/iso,/nod
;;plot, [-1.5,2.5],[-2.5,2.5],/iso,/nod
;;plot, [1.2,2.2], [-0.4,0.4],/nod,/xst,/yst,/iso


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









box=p.(1)
minval=MIN(box.f2d)
maxval=MAX(box.f2d)

for i=1,ntags-1 do begin
   box=p.(i)
   minval=MIN([minval,MIN(box.f2d)])
   maxval=MAX([maxval,MAX(box.f2d)])
endfor

print, "Data range:", minval, maxval


for i=1,ntags-1 do begin

   ;;print, tags[i]
   box=p.(i)
   ;;plot, box.x2d, box.y2d, /iso, psym=3,/xst,/yst
   Polyfill_subgrid, box, minval=minval, maxval=maxval
   ;;STOP
endfor


;
;
;
if keyword_set(DEBUG) then STOP
end
