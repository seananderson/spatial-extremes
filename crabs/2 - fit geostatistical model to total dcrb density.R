library(dplyr)
library(TMB)
#devtools::install_github("nwfsc-assess/geostatistical_delta-GLMM")
library(SpatialDeltaGLMM)

rm(list=ls())

# Calculate spatial information for SPDE mesh, strata areas, and AR1 process
# This is where all runs will be located
DateFile = paste(getwd(),"/crabs/tmb_density",sep="")
# clean out directory
file.remove(paste0(DateFile,"/", dir("crabs/tmb_density")))
if(!exists(DateFile)) dir.create(DateFile)

# All these settings are from WCGTS
Version = "geo_index_v4a"
Method = c("Grid", "Mesh")[2]
grid_size_km = 10
n_x = 300 # Number of stations
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=1, "Epsilon2"=1) # 1=Presence-absence; 2=Density given presence; #Epsilon=Spatio-temporal; #Omega=Spatial
RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0) # Structure for beta or epsilon over time: 0=None (default); 1=WhiteNoise; 2=RandomWalk; 3=Constant
VesselConfig = c("Vessel"=0, "VesselYear"=0)
ObsModel = 2  # 0=normal (log-link); 1=lognormal; 2=gamma; 4=ZANB; 5=ZINB; 11=lognormal-mixture; 12=gamma-mixture
Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )     # Samples: Do K-means on trawl locs; Domain: Do K-means on extrapolation grid

Region = "California_current"
strata.limits <- data.frame(
  'STRATA' = c("Coastwide","CA","OR","WA"),
  'north_border' = c(49.0, 42.0, 46.0, 49.0),
  'south_border' = c(32.0, 32.0, 42.0, 46.0),
  'shallow_border' = c(55, 55, 55, 55),
  'deep_border' = c(1280, 1280, 1280, 1280)
)

dat = readRDS("crabs/data/joined_crab.rds")
dat$cpue = dat$cpue_kgperha * dat$area_swept_ha

Data_Geostat = data.frame( "Catch_KG"=dat[,'cpue'],
  "Year"=dat[,'year'],
  "Vessel"=1,
  "AreaSwept_km2"=dat[,'area_swept_ha'],
  "Lat"=dat[,'lat'],
  "Lon"=dat[,'lon'],
  "Pass"=1,
  "Depth" = dat[,"bottom_depth"]/1000,
  "Depth2" = (dat[,"bottom_depth"]/1000)^2)
Data_Geostat = na.omit( Data_Geostat )

Extrapolation_List = Prepare_Extrapolation_Data_Fn( Region=Region, strata.limits=strata.limits )

Spatial_List = Spatial_Information_Fn( grid_size_km=grid_size_km, n_x=n_x,
  Method=Method, Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'],
  Extrapolation_List=Extrapolation_List, randomseed=Kmeans_Config[["randomseed"]],
  nstart=Kmeans_Config[["nstart"]], iter.max=Kmeans_Config[["iter.max"]], DirPath=DateFile )

Data_Geostat = cbind( Data_Geostat, Spatial_List$loc_UTM, "knot_i"=Spatial_List$knot_i )

# Make TMB data list -- EW added depth as quadratic term in catchability
TmbData = Data_Fn("Version"=Version, "FieldConfig"=FieldConfig, "RhoConfig"=RhoConfig,
  "ObsModel"=ObsModel, "b_i"=Data_Geostat[,'Catch_KG'], "a_i"=Data_Geostat[,'AreaSwept_km2'],
  "Q_ik" = as.matrix(Data_Geostat[,c('Depth','Depth2')]),
  "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1, "s_i"=Data_Geostat[,'knot_i']-1,
  "t_i"=Data_Geostat[,'Year'], "a_xl"=Spatial_List$a_xl, "MeshList"=Spatial_List$MeshList,
  "GridList"=Spatial_List$GridList, "Method"=Spatial_List$Method )

# Make TMB object
TmbList = Build_TMB_Fn("TmbData"=TmbData, "RunDir"=DateFile, "Version"=Version,
  "RhoConfig"=RhoConfig, "VesselConfig"=VesselConfig, "loc_x"=Spatial_List$loc_x)
Obj = TmbList[["Obj"]]

# Run model - EW changed getsd to FALSE
Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]],
  getsd=FALSE, savedir=DateFile, bias.correct=TRUE )
Report = Obj$report()

# Save stuff
Save = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData)
save(Save, file=paste0(DateFile,"Save.RData"))

# PlotAniso_Fn( FileName=paste0(DateFile,"Aniso.png"), Report=Report, TmbData=TmbData )
# Plot surface
#Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
#Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))
#Dim = c( "Nrow"=ceiling(sqrt(length(Years2Include))), "Ncol"=ceiling(length(Years2Include)/ceiling(sqrt(length(Years2Include)))) )
#par( mfrow=Dim )
MapDetails_List = MapDetails_Fn( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )
#PlotResultsOnMap_Fn(plot_set=3, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=paste0(DateFile,"Field_"), Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], mfrow=Dim, mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), Cex=MapDetails_List[["Cex"]], cex=1.8)

# x2i here is the nn.idx index in Spatial_List$PolygonList$NN_Extrap, index of nearest neighbors
# from call to RANN::nn2()
PlotDF = MapDetails_List[["PlotDF"]]
Mat = log(Report$D_xt) # log density at knot locations
incl = which( PlotDF[,'Include']>0 )
Mat = Mat[PlotDF[incl,'x2i'],,drop=FALSE] # expands to nearest neighbors
pred_lon = PlotDF[incl,'Lon']
pred_lat = PlotDF[incl,'Lat']

df = data.frame("lat" = rep(pred_lat, dim(Mat)[2]),
  "lon" = rep(pred_lon, dim(Mat)[2]),
  "logdens" = c(Mat),
  "year" = 2002 + sort(rep(1:13,dim(Mat)[1])))

library(ggplot2)
png("crabs.png", width=700,height=700)
ggplot(data=df, aes(x=lon, y=lat, colour = logdens)) +
  geom_point(size=0.2,alpha=0.2) + facet_wrap(~ year)
dev.off()

Mat.diff = Mat
for(i in 2:dim(Mat.diff)[2]) {
  Mat.diff[,i] = Mat[,i] - Mat[,1]
}
df.diff = data.frame("lat" = rep(pred_lat, dim(Mat)[2]),
  "lon" = rep(pred_lon, dim(Mat)[2]),
  "logdens" = c(Mat.diff),
  "year" = 2002 + sort(rep(1:13,dim(Mat)[1])))

png("crabs/crabs_diff.png", width=700,height=700)
ggplot(data=df.diff[df.diff$year>2003,], aes(x=lon, y=lat, colour = logdens)) +
  geom_point(size=0.2,alpha=0.2) + facet_wrap(~ year)
dev.off()
