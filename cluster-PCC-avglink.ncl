model_name = (/\
"CESM2","IPSL-CM6A-LR","COSMOS","EC-Earth3.3",\
"CESM1.2","IPSLCM5A","MIROC4m","IPSLCM5A2",\
"HadCM3","GISS-E2-1-G","CCSM4","CCSM4-Utr",\
"CCSM4-UoT","NorESM-L","MRI-CGCM2.3", "NorESM1-F","HadGEM3"\
/)
;model_name = (/\
;"CESM2","IPSL-CM6A-LR","EC-Earth3.3",\
;"CESM1.2","IPSLCM5A","MIROC4m","IPSLCM5A2",\
;"HadCM3","GISS-E2-1-G","CCSM4","CCSM4-Utr",\
;"CCSM4-UoT","NorESM-L","MRI-CGCM2.3", "NorESM1-F","HadGEM3"\
;/)
var="pr"
num_m = dimsizes(model_name)
;;;;read data
f         = addfile("cluster_"+var+".nc","r")
distance  = f->distance
distance@_FillValue = 0.0
do i = 0,dimsizes(distance(0,:))-1
 ;print("should be quite small :"+distance(i,i))
 distance(i,i) = 0.0
end do

grp = new((/num_m,num_m/),integer,"999")
grp(0,:) = ispan(0,num_m-1,1)
grp_trace = new((/num_m,num_m,num_m/),integer) ;;(|step,|group matrix)
plt_trace = new((/num_m,3/),integer)

step     = 1
dis      = distance
dis_step = (/0.0/)      ;smallest distance on every step

;;;loop to cluster models
do while(num(.not.ismissing(grp(0,:))).gt.1)
;print(grp(0,:))
  ;;step x find the nearesr 2 groups
  dis1D    = ndtooned(dis)
  dsizes   = dimsizes(dis)
  indices  = ind_resolve(minind(dis1D),dsizes)
  dis_step:= array_append_record(dis_step,min(dis1D),0)
  x        = indices(0,0)
  y        = indices(0,1)
  qx       = num(.not.ismissing(grp(:,x)))        ;; calculate how many members in group x
  qy       = num(.not.ismissing(grp(:,y)))        ;; calculate how many members in group y1
  ;print("this is step: "+step+" and joined member "+x+" member "+y)
  do i = 0,qy-1
  grp(qx+i,x) = grp(i,y)                         ;;append y to group(step)
  end do
  grp(:,y) = grp@_FillValue
  ;;;;;;;;;;;;calculate new distance
  dis           = distance             ;;;reset distance between groups during every step
  if(step.eq.1) then
    rmy = (/y/)
  else
  rmy          := array_append_record(rmy,y,0)
  end if
  dis(rmy,:)    = dis@_FillValue          ;; clear the distance value on grouped model
  dis(:,rmy)    = dis@_FillValue
  do i = 0,num_m-1
   do j = 0,num_m-1
    if(.not.ismissing(dis(i,j))) then
      g2g_dis   = new((/num_m*7/),"float")                 ;; calculate group to group distance,avg of all the distance between members
      mem_i     = grp(:,i)                     ;; indentify two groups member
      mem_j     = grp(:,j)
      do m = 0,num(.not.ismissing(mem_i))-1
       do n = 0,num(.not.ismissing(mem_j))-1
       g2g_dis(m*(num(.not.ismissing(mem_j))+1)+n) = distance(mem_i(m),mem_j(n))
       end do
      end do
      dis(i,j)  = avg(g2g_dis)
    else
      dis(i,j) = dis@_FillValue
    end if
    end do
  end do
  grp_trace(step-1,:,:) = transpose(grp)
  plt_trace(step-1,0)   = step-1
  plt_trace(step-1,1)   = x
  plt_trace(step-1,2)   = y
  step = step+1
end do

do i=0,num_m-1
if(.not.ismissing(grp_trace(step-2,i,0))) then
x_order = grp_trace(step-2,i,:)
end if
end do
xlabel = new((/num_m/),string)
do i = 0,num_m-1
xlabel(i) = model_name(x_order(i))
end do

plt_trace2 = new((/num_m,3/),integer)
plt_trace2(:,0) = plt_trace(:,0)
do i = 0,num_m-1
do j = 0,num_m-1
if(plt_trace(i,1).eq.x_order(j)) then
  plt_trace2(i,1) = j
end if
if(plt_trace(i,2).eq.x_order(j)) then
  plt_trace2(i,2) = j
end if
end do
end do

;;;;start to draw
wks = gsn_open_wks("png","cluster_"+var)

res = True
res@gsnDraw        = False
res@gsnFrame       = False
res@tiMainString   = var
res@vpYF = 0.9
res@tiYAxisString   = "distance"
res@tmXBLabelAngleF  = 90
res@tmXBLabelJust    = "CenterRight"
res@tmXBMode       = "Explicit"         ; explicit labels
res@tmXBValues     = ispan(0,num_m-1,1)
res@tmXTOn         = False
res@tmYROn         = False
res@trXMaxF        = num_m+0.5
res@trXMinF        = -1
res@tmXBLabels     = xlabel
res@trYMaxF        = 1.4
res@trYMinF        = 0.0
y :=new((/num_m/),integer)
plot               = gsn_csm_xy (wks,res@tmXBValues,y,res) ; create plot frame ork

; polygon resources
  pgres                  = True
  pgres@gsLineThicknessF = 2.0   ; make lines thicker
  dumb1  = new(step-1,graphic)
  dumb2  = new(step-1,graphic)
  dumb3  = new(step-1,graphic)
  h_length = new((/step-1/),"float")
  do pstep = 0,step-2
    x_sta = int2flt(plt_trace2(pstep,1))
    x_end = int2flt(plt_trace2(pstep,2))
    ;a=num(.not.ismissing(grp_trace(pstep,plt_trace(pstep,1),:))) ;calculate how many members under a group
  if(pstep.gt.1) then
    do i = 0,pstep-1
      if(any(plt_trace2(i,1:2).eq.plt_trace2(pstep,1))) then
      ;x_sta = x_sta+0.5
      x_sta = x_sta+0.5*h_length(i)
      end if
    end do
    do i = 0,pstep-1
      if(any(plt_trace2(i,1:2).eq.plt_trace2(pstep,2))) then
      ;x_end = x_end+0.5
      x_end = x_end+0.5*h_length(i)
      end if
    end do
  end if
  h_length(pstep)=x_end-x_sta
  y_sta1 = 0.0
  y_sta2 = 0.0
    do i = 0,pstep-1
      if(any(plt_trace2(i,1:2).eq.plt_trace2(pstep,1))) then
      y_sta1 = dis_step(i+1)
      end if
    end do
    do i = 0,pstep-1
      if(any(plt_trace2(i,1:2).eq.plt_trace2(pstep,2))) then
      y_sta2 = dis_step(i+1)
      end if
    end do
  y_end = dis_step(pstep+1)
  dumb1(pstep) = gsn_add_polyline(wks,plot,(/x_sta,x_end/),(/y_end,y_end/),pgres)
  dumb2(pstep) = gsn_add_polyline(wks,plot,(/x_sta,x_sta/),(/y_sta1,y_end/),pgres)
  dumb3(pstep) = gsn_add_polyline(wks,plot,(/x_end,x_end/),(/y_sta2,y_end/),pgres)
  end do

;;-----draw step fig
 ;wks2   = gsn_open_wks ("png","xy")                  ; send graphics to PNG file
 res2                   = True                       ; plot mods desired
 res2@gsnDraw        = False
 res2@gsnFrame       = False
 res2@vpWidthF        = .20          ; set width of second plot
 res2@xyMarkLineModes   = "Markers"                ; choose which have markers
 res2@xyMarkers         =  16                     ; choose type of marker  
 res2@tmXBLabels     = ispan(1,step-1,1)
 res2@trYMaxF        = res@trYMaxF
 res2@trYMinF        = res@trYMinF
 res2@tiYAxisString   = "distance"
 res2@tiXAxisString   = "step"
 plot2  = gsn_csm_xy (wks,res2@tmXBLabels,dis_step(1::),res2) ; create plot

 newplot = gsn_attach_plots(plot,(/plot2/),res,res2)

draw(plot)
frame(wks)
