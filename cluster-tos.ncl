
var_ind = 0
sea_ind = 0
var_name = (/"tos"/)
;model_name = (/\
;"CESM2","IPSL-CM6A-LR","COSMOS","EC-Earth3.3",\
;"CESM1.2","IPSLCM5A","MIROC4m","IPSLCM5A2",\
;"HadCM3","GISS-E2-1-G","CCSM4","CCSM4-Utr",\
;"CCSM4-UoT","NorESM-L","MRI-CGCM2.3", "NorESM1-F","HadGEM3"\
;/)
model_name = (/\
"CESM2","IPSL-CM6A-LR","EC-Earth3.3",\
"CESM1.2","IPSLCM5A","MIROC4m","IPSLCM5A2",\
"HadCM3","GISS-E2-1-G","CCSM4","CCSM4-Utr",\
"CCSM4-UoT","NorESM-L","MRI-CGCM2.3", "NorESM1-F","HadGEM3"\
/)
;;-----set latlon box
  latt = 15.
  latb = -15.
  lonl = 80.
  lonr = 160.
;;------set the final var
  fi  = ("/export/silurian/array-01/bo20541/processed_data/tos/tos_regrid_50y_12month/"+var_name(var_ind)+"_"+model_name(10)+"_Omon_Eoi400_gr_01-50.nc")
  fii = addfile(fi,"r")
  dem = fii->tos
  lat = fii->lat
  lon = fii->lon
  demo= dem(0,{latb:latt},{lonl:lonr})
  var = new((/dimsizes(model_name),dimsizes(demo(:,0)),dimsizes(demo(0,:))/),typeof(demo),demo@_FillValue)
  var!0 = "model_name"
  var!1 = "lat"
  var!2 = "lon"
  var@model_name = model_name
  var&lat = lat({latb:latt})
  var&lon = lon({lonl:lonr})

;;----read data from mdoels
do mod_i = 0,dimsizes(model_name)-1
 test_file  =("/home/bridge/bo20541/work/processed_data/"+var_name(var_ind)+"/tos_regrid_50y_12month/"+var_name(var_ind)+"_"+model_name(mod_i)+"_Omon_Eoi400_gr_01-50.nc")
if(fileexists(test_file)) then
  ;;--------------------------------------------------------------------
  file_path1  = ("/export/silurian/array-01/bo20541/processed_data/tos/tos_regrid_50y_12month/"+var_name(var_ind)+"_"+model_name(mod_i)+"_Omon_Eoi400_gr_01-50.nc")
  file_path0  = ("/export/silurian/array-01/bo20541/processed_data/tos/tos_regrid_50y_12month/"+var_name(var_ind)+"_"+model_name(mod_i)+"_Omon_E280_gr_01-50.nc")
  f1          = addfile(file_path1,"r")
  f0          = addfile(file_path0,"r")
  var1        = f1->$(var_name(var_ind))$
  var0        = f0->$(var_name(var_ind))$
  if(model_name(mod_i).eq."HadGEM3") then
   var1@_FillValue = getVarFillValue(var1)
   var1            = where(var1.gt.-200.and.var1.lt.200,var1,var1@_FillValue)
   var0            = where(var0.gt.-200.and.var0.lt.200,var0,var1@_FillValue)
  end if
  vard        = var1 - var0
  copy_VarCoords(var1,vard)
  dataANN     = dim_avg_n(vard(:,{latb:latt},{lonl:lonr}),0)
  var(mod_i,:,:) = dataANN/stddev(dataANN)
  delete([/var1,var0,vard,dataANN/])
else
  var(mod_i,:,:) = demo@_FillValue
end if
end do

 gw   = latRegWgt(lat({latb:latt}), "float", 0)
 r    = new((/dimsizes(model_name),dimsizes(model_name)/),"float")
;;;for model x
do i = 0,dimsizes(model_name)-1
  do j = 0,dimsizes(model_name)-1
           r(i,j) = pattern_cor(var(i,:,:), var(j,:,:),  gw, 0)      ; gaussian weighting, centered
  end do
end do
r=1.0-r

fo=addfile("pliomip/plot_script/cluster_tos.nc","c")
fo=addfile("pliomip/plot_script/cluster_tos.nc","w")
fo->distance=r
