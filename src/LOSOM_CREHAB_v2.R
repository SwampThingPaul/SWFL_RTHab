## 
## LOSOM Water Quality Subteam
## Red tide - HAB model
##
## Code was compiled by Paul Julian
## contact info: paul.julian@floridadep.gov

## BAD 
## https://www.tidyverse.org/articles/2017/12/workflow-vs-script/
## Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

## Libraries
#devtools::install_github("SwampThingPaul/AnalystHelper")
library(AnalystHelper);
library(plyr)
library(reshape)
library(zoo)

# GIS libraries 
library(rgdal)
library(rgeos)
library(raster)
library(tmap)

library(classInt)
library(mgcv)
library(lubridate)
library(gratia)
library(ggplot2)

library(flextable)
library(magrittr)

## Paths
wd="C:/Julian_LaCie/_Github/LOSOM_RTHab"

paths=paste0(wd,c("/Plots/","/Export/","/Data/","/GIS","/src/"))
# Folder.Maker(paths);#One and done. Creates folders in working directory.
plot.path=paths[1]
export.path=paths[2]
data.path=paths[3]
GIS.path=paths[4]

GIS.path="C:/Julian_LaCie/_GISData"

## Functions
notidy_glance_gam<-function(model,...){
  data.frame(
    df=sum(model$edf),
    df.residual=stats::df.residual(model),
    logLik=as.numeric(stats::logLik(model)),
    AIC = stats::AIC(model),
    BIC = stats::BIC(model),
    adj.r.squared=summary(model)$r.sq,
    deviance=summary(model)$dev.expl,
    nobs = stats::nobs(model),
    method=as.character(summary(model)$method),
    sp.crit=as.numeric(summary(model)$sp.criterion),
    scale.est=summary(model)$scale
  )
}

notidy_tidy_gam<-function(model,dig.num=2,...){
  ptab <- data.frame(summary(model)$p.table)
  ptab$term<-rownames(ptab)
  rownames(ptab)=NULL
  ptab$Component="A. parametric coefficients"
  ptab<-ptab[,c(6,5,1:4)]
  colnames(ptab) <- c("Component","Term", "Estimate", "Std.Error", "t.value", "p.value")
  ptab$p.value=with(ptab,ifelse(p.value<0.01,"<0.01",round(p.value,2)))
  ptab[,3:5]=format(round(ptab[,3:5],dig.num),nsmall=dig.num)
  ptab
  
  stab= data.frame(summary(model)$s.table)
  stab$term<-rownames(stab)
  rownames(stab)=NULL
  stab$Component="B. smooth terms"
  stab<-stab[,c(6,5,1:4)]
  colnames(stab) <- c("Component","Term", "edf", "Ref. df", "F.value", "p.value")
  stab$p.value=with(stab,ifelse(p.value<0.01,"<0.01",round(p.value,2)))
  stab[,3:5]=format(round(stab[,3:5],dig.num),nsmall=dig.num)
  stab
  
  ptab.cnames = c("Component","Term", "Estimate", "Std Error", "t-value", "p-value")
  stab.cnames = c("Component","Term", "edf", "Ref. df", "F-value", "p-value")
  
  colnames(ptab) = c("A", "B", "C", "D")
  if (ncol(stab) != 0) {
    colnames(stab) = colnames(ptab)
  }
  tab = rbind(ptab, stab)
  colnames(tab) = ptab.cnames
  
  tab2 = rbind(c(ptab.cnames), tab[1:nrow(ptab), ])
  if (nrow(stab) > 0) {
    tab2 = rbind(tab2, c(stab.cnames), tab[(nrow(ptab) + 1):nrow(tab), ])
  }
  
  tab2
}

notidy_as_flextable_gam<-function(x,data_t=NULL,data_g=NULL,dig.num=2,r2dig=2,...){
  # needs flextable
  # magrittr
  if(sum(class(x)%in%c("gam"))==1&is.null(data_t)&is.null(data_g)){
    data_t <- notidy_tidy_gam(x)
    data_g <- notidy_glance_gam(x)
  }
  
  std_border=officer::fp_border(color = "black", style = "solid", width = 2)
  data.frame(data_t)%>%
    flextable()%>%
    delete_part(part="header")%>%
    hline(i=which(data_t=="Component"),border=std_border)%>%
    hline(i=which(data_t=="Component")[2]-1,border=std_border)%>%
    bold(i=which(data_t=="Component"))%>%
    align(j=1,part="all")%>%
    hline_top(border=std_border)%>%
    hline_bottom(border=std_border)%>%
    merge_v(j=1)%>%valign(j=1,valign="top")%>%fix_border_issues()%>%
    autofit(part = c("header", "body"))%>%
    add_footer_lines(values = c(
      sprintf("Adjusted R-squared: %s, Deviance explained %s", formatC(data_g$adj.r.squared,digits = r2dig,format="f"), formatC(data_g$deviance,digits = r2dig,format="f")),
      paste0(data_g$method,": ",format(round(data_g$sp.crit,dig.num),dig.num),", Scale est.: ",format(round(data_g$scale.est,dig.num),dig.num),", N: ",data_g$nobs)
    ))
}

nad83.pro=CRS("+init=epsg:4269")
utm17=CRS("+init=epsg:26917")

tmap_mode("view")
# -------------------------------------------------------------------------
# GIS ---------------------------------------------------------------------
shore=spTransform(readOGR(paste0(GIS.path,"/FWC"),"FWC_Shoreline"),wkt(utm17))
shore2=gSimplify(shore,250)

wmd.struct=spTransform(readOGR(paste0(GIS.path,"/AHED_release/AHED_20171102.gdb"),"STRUCTURE"),wkt(utm17))
canals=spTransform(readOGR(paste0(GIS.path,"/SFER_GIS_Geodatabase.gdb"),"SFWMD_Canals"),wkt(utm17))
roads.all=spTransform(readOGR(paste0(GIS.path,"/FDOT"),"FDOT_Roads"),wkt(utm17))
lakes=spTransform(readOGR(paste0(GIS.path,"/NHD"),"NHD100_Waterbody"),wkt(utm17))
wetland=subset(lakes,FTYPE%in%c("466"))

## AOI (general area from Medina et al.)
bbox=raster::extent(-82.5,-81.65,26,27)
bbox.poly=as(bbox,"SpatialPolygons")
proj4string(bbox.poly)=nad83.pro
bbox.poly=spTransform(bbox.poly,utm17)
bbox.poly=SpatialPolygonsDataFrame(bbox.poly,data.frame(ID=1))

shore.clip=raster::intersect(gSimplify(shore,500),bbox.poly)
shore.clip.buf=gBuffer(shore.clip,width=16*1000)

AOI=raster::intersect(gBuffer(bbox.poly,width=0.5),shore.clip.buf)
plot(AOI)
##

example=c("01-APR-02 05.45.00.000000000 PM", "01-APR-02 12.00.00.000000000 AM", 
          "01-APR-02 12.00.00.000000000 AM", "01-APR-02 12.00.00.000000000 AM", 
          "01-APR-02 12.00.00.000000000 AM", "01-APR-02 12.00.00.000000000 AM")
spl=strsplit(example,split="-")
date_split=data.frame(SAMPLE_DATE=example,
                      day=sapply(spl,"[",1),month=sapply(spl,"[",2),yr_time=sapply(spl,"[",3))
date_split2=strsplit(date_split$yr_time,split=" ")
date_split=cbind(date_split[,1:3],data.frame(yr=sapply(date_split2,"[",1),time=sapply(date_split2,"[",2),AMPM=sapply(date_split2,"[",3)))
date_split$yr=with(date_split,ifelse(yr%in%seq(53,99,1),paste0("19",yr),paste0("20",yr)))
date_split$SAMPLE_DATE2=with(date_split,paste(paste(day,month,yr,sep="-"),time,AMPM))
date_split$datetime=date.fun(date_split$SAMPLE_DATE2,form="%d-%b-%Y %I.%M.%OS %p")
date_split[,c("SAMPLE_DATE","datetime")]

# HAB data ----------------------------------------------------------------
# hab.dat=read.csv(paste0(data.path,"habsos_20200310.csv"))
hab.dat=read.csv(paste0(data.path,"habsos_20210413.csv")); #new data file

spl=strsplit(hab.dat$SAMPLE_DATE,split="-")
date_split=data.frame(SAMPLE_DATE=hab.dat$SAMPLE_DATE,
                 day=sapply(spl,"[",1),month=sapply(spl,"[",2),yr_time=sapply(spl,"[",3))
date_split2=strsplit(date_split$yr_time,split=" ")
date_split=cbind(date_split[,1:3],data.frame(yr=sapply(date_split2,"[",1),time=sapply(date_split2,"[",2),AMPM=sapply(date_split2,"[",3)))
date_split$yr=with(date_split,ifelse(yr%in%seq(53,99,1),paste0("19",yr),paste0("20",yr)))
date_split$SAMPLE_DATE2=with(date_split,paste(paste(day,month,yr,sep="-"),time,AMPM))
date_split$datetime=date.fun(date_split$SAMPLE_DATE2,form="%d-%b-%Y %I.%M.%OS %p")

# hab.dat$datetime=date.fun(hab.dat$SAMPLE_DATE,form="%d-%b-%y %I.%M.%OS %p")
hab.dat=merge(hab.dat,date_split[,c("SAMPLE_DATE","datetime")],"SAMPLE_DATE")
hab.dat$date=date.fun(hab.dat$datetime)
hab.dat$biweekperiod=as.numeric(format(hab.dat$date,"%j"))%/%14L+1L
hab.dat$month=as.numeric(format(hab.dat$date,'%m'))
hab.dat$CY=as.numeric(format(hab.dat$date,'%Y'))

subset(hab.dat,date==date.fun("2068-05-08"))

hab.dat.shp=spTransform(SpatialPointsDataFrame(hab.dat[,c("LONGITUDE","LATITUDE")],data=hab.dat,proj4string=nad83.pro),utm17)
hab.dat.shp.AOI=hab.dat.shp[AOI,]

# png(filename=paste0(plot.path,"HABSOS_map.png"),width=5,height=6,units="in",res=200,type="windows",bg="white")
# par(family="serif",mar=c(0.5,0.5,0.5,0.5),oma=c(1,1,1,1));
# layout(matrix(1:2,1,2,byrow=F),widths=c(1,0.5))

# bbox.lims=bbox(gBuffer(AOI,width=1000))
# plot(shore2,ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],col="cornsilk",border="grey",bg="lightblue",lwd=0.2,xpd=F)
# plot(lakes,add=T,border=NA,col=adjustcolor("skyblue",0.5))
# plot(wetland,add=T,border=NA,col=adjustcolor("palegreen3",0.5))
# plot(canals,add=T,col="dodgerblue1",lwd=2)
# plot(canals,add=T,col="lightblue",lwd=1)
# plot(roads.all,add=T,col="grey",lwd=0.75,lty=1)
# plot(wmd.struct,pch=19,col="black",add=T)
# plot(AOI,add=T,border="red",lty=2,lwd=2)
# plot(subset(hab.dat.shp,CELLCOUNT>0),add=T,pch=21,bg=adjustcolor("white",0.25),col=adjustcolor("grey",0.25))
# plot(subset(hab.dat.shp.AOI,CELLCOUNT>0),add=T,pch=21,bg=adjustcolor("green",0.25),col=adjustcolor("black",0.25))
# mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=0.75,seg.len=4);box(lwd=1)
# 
# # plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA);
# leg.text=c("AOI", "HABSOS Data - All", "HABSOS Data - Study Area")
# legend("bottomleft",legend=leg.text,
#        pch=c(NA,21,21),
#        lty=c(2,NA,NA),lwd=c(2,0.5,0.5),
#        col=c("red",adjustcolor(c("grey","black"),0.25)),pt.bg=c(NA,adjustcolor(c("white","green"),0.25)),
#        pt.cex=2,ncol=1,cex=0.8,bty="n",y.intersp=1.1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
# dev.off()

hab.dat2=hab.dat.shp.AOI@data

# fill.days=seq(date.fun("1953-05-01"),date.fun("2020-05-01"),"1 days")
# fill.weeks=as.numeric(format(fill.days,"%j"))%/%14L+1L
# hab.dat3=merge(hab.dat2,data.frame(date=fill.days,CY=as.numeric(format(fill.days,'%Y')),biweekperiod=fill.weeks,fill=1),c("date","CY","biweekperiod"),all.y=T)

hab.biweek=ddply(subset(hab.dat2,SAMPLE_DEPTH<=5),c("CY","biweekperiod"),summarise,
                 min.date=min(date,na.rm=T),
                 max.date=max(date,na.rm=T),
                 mean.val=mean(CELLCOUNT,na.rm=T),
                 med.val=median(CELLCOUNT,na.rm=T),
                 min.val=min(CELLCOUNT,na.rm=T),
                 max.val=max(CELLCOUNT,na.rm=T))
hab.biweek$log10.mean=with(hab.biweek,ifelse(mean.val==0,0,log10(mean.val)))

ylim.val=c(0,7);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=date.fun(c("1953-05-01","2021-02-24"));xmaj=seq(xlim.val[1],xlim.val[2],"10 years");xmin=seq(xlim.val[1],xlim.val[2],"6 months")
# png(filename=paste0(plot.path,"HAB_TS_all.png"),width=7,height=2.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,4,0.25,0.25),oma=c(2,2,0.75,0.5));

plot(log10.mean~min.date,hab.biweek,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(hab.biweek,segments(min.date,0,min.date,log10.mean,col=adjustcolor("black",0.5)))
with(hab.biweek,points(min.date,log10.mean,pch=19,cex=0.5,col=adjustcolor("black",0.5)))
axis_fun(1,xmaj,xmin,format(xmaj,"%m-%Y"),line=-0.5)
axis_fun(2,ymaj,ymaj,format(c(0,10^(ymaj[2:8])),scientific = F));box(lwd=1)
mtext(side=2,line=4,expression(paste(italic("K. brevis")," (cells L"^"-1",")")))
mtext(side=1,line=1.75,"Date (Month-Year)")
abline(v=date.fun(c("2012-12-01","2018-02-28")),lty=3,col="red")

dev.off()

ylim.val=c(0,7);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=date.fun(c("2012-12-01","2021-02-24"));xmaj=seq(xlim.val[1],xlim.val[2],"10 years");xmin=seq(xlim.val[1],xlim.val[2],"6 months")
par(family="serif",mar=c(1,4,0.25,0.25),oma=c(2,2,0.75,0.5));

plot(log10.mean~min.date,hab.biweek,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(hab.biweek,segments(min.date,0,min.date,log10.mean,col=adjustcolor("black",0.5)))
with(hab.biweek,points(min.date,log10.mean,pch=19,cex=0.5,col=adjustcolor("black",0.5)))
axis_fun(1,xmaj,xmin,format(xmaj,"%m-%Y"),line=-0.5)
axis_fun(2,ymaj,ymaj,format(c(0,10^(ymaj[2:8])),scientific = F));box(lwd=1)
mtext(side=2,line=4,expression(paste(italic("K. brevis")," (cells L"^"-1",")")))
mtext(side=1,line=1.75,"Date (Month-Year)")
abline(v=date.fun(c("2012-12-01","2018-02-28")),lty=3,col="red")



# Discharge ---------------------------------------------------------------
dates=date.fun(c("1995-05-01","2021-12-31"))

q.dbkeys=data.frame(SITE=c("S79","S78","S77","S65E","S65EX1"),DBKEY=c("00865","DJ236","15635","91656","AL760"))
q.dat=data.frame()
for(i in 1:nrow(q.dbkeys)){
  tmp=DBHYDRO_daily(dates[1],dates[2],q.dbkeys$DBKEY[i])
  tmp$DBKEY=as.character(q.dbkeys$DBKEY[i])
  q.dat=rbind(q.dat,tmp)
  print(i)
}
q.dat=merge(q.dat,q.dbkeys,"DBKEY")
q.dat$Date.EST=date.fun(q.dat$Date)
q.dat$WY=WY(q.dat$Date.EST)

q.dat.xtab=cast(q.dat,Date.EST~SITE,value="Data.Value",mean)

q.dat.xtab$month=as.numeric(format(q.dat.xtab$Date.EST,"%m"))
q.dat.xtab$CY=format(q.dat.xtab$Date.EST,"%Y")
q.dat.xtab$monCY=with(q.dat.xtab,date.fun(paste(CY,month,"01",sep="-")))
q.dat.xtab$WY=WY(q.dat.xtab$Date.EST)

#Positive discharge only
q.dat.xtab$S77=with(q.dat.xtab,ifelse(S77<0,0,S77))
q.dat.xtab$S79=with(q.dat.xtab,ifelse(S79<0,0,S79))
q.dat.xtab$C43=with(q.dat.xtab,ifelse(S79<S77,0,S79-S77))
q.dat.xtab$S65E_T=rowSums(q.dat.xtab[,c("S65E","S65EX1")],na.rm=T)
q.dat.xtab$biweekperiod=as.numeric(format(q.dat.xtab$Date.EST,"%j"))%/%14L+1L
q.dat.xtab$Date.EST.c=as.character(q.dat.xtab$Date.EST)
head(q.dat.xtab)

vars=c("Date.EST","month","CY","biweekperiod","S65E_T","S77","S79","C43")
q.dat.xtab.melt=melt(data.frame(q.dat.xtab[,vars]),id.vars=vars[1:4])

q.dat.xtab2=cast(q.dat.xtab.melt,CY+biweekperiod~variable,value="value",fun.aggregate = function(x)sum(x,na.rm=T))
q.dat.xtab.mon=cast(q.dat.xtab.melt,CY+month~variable,value="value",fun.aggregate = function(x)sum(x,na.rm=T))

plot(S79~Date.EST,q.dat.xtab)
plot(C43~Date.EST,q.dat.xtab)

# WQ ----------------------------------------------------------------------
params=data.frame(Test.Number=c(18,21,80,20,25,23,61,179,7,16),
                  param=c("NOx","TKN","TN","NH4","TP","SRP","Chla","Chla","Temp","TSS"))
params=subset(params,param%in%c("NOx","TKN","TN","TP"))
wq.sites="S79" #c("S79","S77","S65E")
wq.dat=DBHYDRO_WQ(dates[1],dates[2],wq.sites,params$Test.Number)
wq.dat=merge(wq.dat,params,"Test.Number")
unique(wq.dat$Collection.Method)
wq.dat=subset(wq.dat,Collection.Method=="G")

# plot(HalfMDL~Date.EST,subset(wq.dat,Station.ID=="S79"& param=="TP"))
# plot(HalfMDL~Date.EST,subset(wq.dat,Station.ID=="S77"& param=="TP"))

wq.dat.xtab=cast(wq.dat,Station.ID+Date.EST~param,value="HalfMDL",mean)
wq.dat.xtab$TN=with(wq.dat.xtab, TN_Combine(NOx,TKN,TN))
wq.dat.xtab$WY=WY(wq.dat.xtab$Date.EST)
wq.dat.xtab$month=as.numeric(format(wq.dat.xtab$Date.EST,"%m"))
wq.dat.xtab$CY=format(wq.dat.xtab$Date.EST,"%Y")
wq.dat.xtab$monCY=with(wq.dat.xtab,date.fun(paste(CY,month,"01",sep="-")))
wq.dat.xtab$biweekperiod=as.numeric(format(wq.dat.xtab$Date.EST,"%j"))%/%14L+1L

vars=c("Station.ID","Date.EST","month","CY","biweekperiod","TN","TP")
wq.dat.xtab.melt=melt(data.frame(wq.dat.xtab[,vars]),id.vars=vars[1:5])
wq.dat.xtab.melt$Station_param=with(wq.dat.xtab.melt,paste(Station.ID,variable,sep="_"))

wq.dat.xtab2=cast(wq.dat.xtab.melt,CY+biweekperiod~Station_param,value="value",fun.aggregate = function(x)mean(x,na.rm=T))

wq.dat.xtab.mon=cast(wq.dat.xtab.melt,CY+month~Station_param,value="value",fun.aggregate = function(x)mean(x,na.rm=T))
plot(TN~Date.EST,subset(wq.dat.xtab,Station.ID=="S79"))

# Lake O Stage ------------------------------------------------------------
lakeO.stg=DBHYDRO_daily(dates[1],dates[2],"00268")
plot(Data.Value~Date,lakeO.stg)

lakeO.stg$Date.EST=date.fun(lakeO.stg$Date)
lakeO.stg$WY=WY(lakeO.stg$Date)
lakeO.stg$CY=as.numeric(format(lakeO.stg$Date,"%Y"))
lakeO.stg$month=as.numeric(format(lakeO.stg$Date,"%m"))
lakeO.stg$biweekperiod=as.numeric(format(lakeO.stg$Date,"%j"))%/%14L+1L

lakeO.stg2=ddply(lakeO.stg,c("CY","biweekperiod"),summarise,mean.stg=mean(Data.Value,na.rm=T))
lakeO.stg.mon=ddply(lakeO.stg,c("CY","month"),summarise,mean.stg=mean(Data.Value,na.rm=T))

lakeO.stg$DoY=hydro.day(lakeO.stg$Date.EST)

#LORS

load("C:/Julian_LaCie/_GitHub/CRE_Conditions/Data/LORS.Rdata")
CurWY=2012
LORS$Year=CurWY-1
LORS$Date2=with(LORS,date.fun(paste(Year,Month,Day,sep="-")))
LORS$DoY=hydro.day(LORS$Date2)
LORS=LORS[order(LORS$DoY),]

HighLakeLab.x=hydro.day(as.POSIXct(paste(CurWY-1,9,1,sep="-")))
WSMLab.x=hydro.day(as.POSIXct(paste(CurWY-1,9,1,sep="-")))
BENLab.x=hydro.day(as.POSIXct(paste(CurWY-1,7,15,sep="-")))
BASELab.x=hydro.day(as.POSIXct(paste(CurWY,4,1,sep="-")))
HIGHLab.x=hydro.day(as.POSIXct(paste(CurWY,3,15,sep="-")))
InterLab.x=hydro.day(as.POSIXct(paste(CurWY,3,1,sep="-")))
LowLab.x=hydro.day(as.POSIXct(paste(CurWY,2,15,sep="-")))

lwd.val=1
xlim.vals=c(date.fun(paste(CurWY-1,05,01,sep="-")),date.fun(paste(CurWY,04,30,sep="-")))
xmaj.lab=seq(xlim.vals[1],xlim.vals[2],by="2 months");xmin.lab=seq(xlim.vals[1],xlim.vals[2],by="1 months")
xlim.vals=hydro.day(xlim.vals)
xmaj=hydro.day(xmaj.lab);xmin=hydro.day(xmin.lab)

ylim.val=c(9,18);by.y=1
ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
# png(filename=paste0(plot.path,"LORS.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(0.5,2,1,1),oma=c(0.5,3,1,1));
layout(matrix(1:2,2,1,byrow=T),heights=c(1,0.3))
plot(High~DoY,LORS,ylim=ylim.val,xlim=xlim.vals,type="n",lwd=2,ylab=NA,xlab=NA,yaxs="i",xaxs="i",xaxt="n",yaxt="n")
abline(h=seq(9,18,1),lwd=1,col="grey",lty=3)
abline(v=xmin,lwd=1,col="grey",lty=3)
with(LORS,lines(High~DoY,lwd=2,col="black"))
with(LORS,lines(Intermediate~DoY,lwd=2,col="black"))
with(LORS,lines(Low~DoY,lwd=2,col="black"))
with(LORS,lines(BaseFlow~DoY,lwd=2,col="black"))
with(LORS,lines(BeneficialUse~DoY,lwd=2,col="black"))
with(LORS,lines(WSM~DoY,lwd=2,col="grey"))
with(LORS,lines(Inter1ft~DoY,lwd=2,lty=5,col="black"))
with(LORS,lines(LowLow~DoY,lwd=2,lty=5,col="grey"))
with(LORS,lines(LowMid~DoY,lwd=2,lty=5,col="grey"))
text(HighLakeLab.x,17.5,"High Lake Management Band",font=2)
text(WSMLab.x,9.5,"Water Shortage Management Band",font=2)
text(BENLab.x,12,"Beneficial Use",font=2)
text(BASELab.x,13,"Base Flow",font=2)
text(HIGHLab.x,17,"High",font=2)
text(InterLab.x,16.25,"Intermediate ",font=2,cex=0.75)
text(LowLab.x,14.5,"Low",font=2,cex=0.75)
axis_fun(1,xmaj,xmin,format(xmaj.lab,"%b"),line=-0.5)
axis_fun(2,ymaj,ymin,ymaj)
box(lwd=lwd.val)
mtext(side=1,"Month",line=2,cex=1.25)
mtext(side=2,"Stage Elevation (Feet, NGVD29)",line=2.5,cex=1.25)
cols=wesanderson::wes_palette("Zissou1",7,"continuous")
Yrs=seq(2012,2018,1)
for(i in 1:7){
  with(subset(lakeO.stg,WY==Yrs[i]),lines(DoY,Data.Value,col=adjustcolor(cols[i],0.5),lwd=3))
}
mtext(side=3,adj=0,"Lake Okeechobee - LORS08")

plot(0:1,0:1,type = 'n', axes = F,xlab=NA, ylab=NA)
legend(0.5,0.5,
       legend=paste0("WY",seq(2012,2018,1)),
       col=adjustcolor(cols,0.5),lty=c(1),lwd=c(4),ncol=4,cex=0.75,bty="n",y.intersp=1,x.intersp=0.5,xpd=NA,xjust=0.5)

dev.off()


# Expanded data -----------------------------------------------------------
expanded.dates=date.fun(c("2012-10-01","2018-10-01"))
# expanded.dates=date.fun(c("2012-10-01","2020-10-01"))
# expanded.dates=date.fun(c("1995-10-01","2010-10-01"))

hab.dat2$log10.mean=with(hab.dat2,ifelse(CELLCOUNT==0,0,log10(CELLCOUNT)))
hab.dat2=hab.dat2[order(hab.dat2$date),]

plot(log10.mean~date,hab.dat2,xlim=expanded.dates,type="l")
abline(v=expanded.dates,col="red",lty=2,lwd=2)


# Salinity ----------------------------------------------------------------
# surf.da=data.frame(depth="surface",SITE=c(rep("S79_T",2),rep("VALI75",2),rep("FORTMYERSM",4),rep("CCORAL",2),rep("MARKH",2),rep("MARKERH",2),rep("SANIB1",2)),param=c(rep(c("WT","SPC"),8)),DBKEY=c("15286","15287","UL031","UL027","PE681","PE685","88285","88289","UO833","PS985","88199","88203","15269","15271","WN366","WN368"))
# bot.da=data.frame(depth="bottom", SITE=c(rep("S79_T",2),rep("VALI75",2),rep("FORTMYERSM",4),rep("CCORAL",2),rep("MARKH",2),rep("MARKERH",2),rep("SANIB2",2)),param=c(rep(c("WT","SPC"),8)),DBKEY=c("15285","15289","UL029","UL025","PE684","PE688","88286","88290","UO831","PS984","88197","88201","15268","15270","WN374","WN376"))
# 
# sal.dat=data.frame()
# for(i in 1:nrow(surf.da)){
#   tmp=DBHYDRO_daily(dates[1],dates[2],surf.da$DBKEY[i])
#   tmp$DBKEY=as.character(surf.da$DBKEY[i])
#   sal.dat=rbind(sal.dat,tmp)
#   print(i)
# }
# sal.dat=merge(sal.dat,surf.da,"DBKEY")
# sal.dat.xtab=data.frame(cast(sal.dat,Date+SITE~param,value="Data.Value",mean))
# sal.dat.xtab$Sal=with(sal.dat.xtab,SalinityCalc(SPC,WT))
# 
# sal.dat.xtab2=data.frame(cast(sal.dat.xtab,Date~SITE,value="Sal",mean))
# sal.dat.xtab2=subset(sal.dat.xtab2,Date>date.fun("2006-10-01"))# VALI75 data starts 2006-10-20
# 
# plot(MARKH~Date,sal.dat.xtab2)
# plot(FORTMYERSM~Date,sal.dat.xtab2)
# 
# sal.dat.xtab2$est.sal.delta=with(sal.dat.xtab2,MARKH-FORTMYERSM)
# sal.dat.xtab2$month=as.numeric(format(sal.dat.xtab2$Date,"%m"))
# sal.dat.xtab2$CY=as.numeric(format(sal.dat.xtab2$Date,"%Y"))
#   
# plot(est.sal.delta~Date,sal.dat.xtab2)
# sal.dat.xtab2.mon=ddply(sal.dat.xtab2,c("CY","month"),summarise,mean.sal.delta=mean(est.sal.delta,na.rm=T),mean.FtM=mean(FORTMYERSM,na.rm=T),mean.MH=mean(MARKH,na.rm=T))
# plot(sal.dat.xtab2.mon$mean.sal.delta)
# plot(sal.dat.xtab2.mon$mean.FtM)

# -------------------------------------------------------------------------
exp.hab=subset(hab.dat2,date>=expanded.dates[1]&date<expanded.dates[2])

# Expanded data - Monthly -------------------------------------------------
exp.hab.month=ddply(subset(exp.hab,SAMPLE_DEPTH<=5),c("CY","month"),summarise,
                min.date=min(date),
                max.date=max(date),
                mean.val=mean(CELLCOUNT,na.rm=T),
                med.val=median(CELLCOUNT,na.rm=T),
                min.val=min(CELLCOUNT,na.rm=T),
                max.val=max(CELLCOUNT,na.rm=T),
                N.val=N.obs(CELLCOUNT))
exp.hab.month$log10.mean=with(exp.hab.month,ifelse(mean.val==0,0,log10(mean.val)))
exp.hab.month$rec=1
exp.hab.month$rec[1]=exp.hab.month$month[1]
exp.hab.month$c.month=cumsum(exp.hab.month$rec)

exp.hab.month=merge(exp.hab.month,lakeO.stg.mon,c("CY","month"),all.x=T)
exp.hab.month=merge(exp.hab.month,data.frame(wq.dat.xtab.mon),c("CY","month"),all.x=T)
exp.hab.month=merge(exp.hab.month,data.frame(q.dat.xtab.mon),c("CY","month"),all.x=T)
exp.hab.month$date.monCY=with(exp.hab.month,date.fun(paste(CY,month,'01',sep="-")))
exp.hab.month$DoY=as.numeric(format(exp.hab.month$date.monCY,"%j"))
exp.hab.month$DY=with(exp.hab.month,lubridate::decimal_date(date.monCY))
# exp.hab.month$DY=with(exp.hab.month,lubridate::decimal_date(date.fun("1965-05-01"))-lubridate::decimal_date(date.monCY))
exp.hab.month$WY=WY(exp.hab.month$min.date,"Fed")
exp.hab.month$DOWY=as.numeric(hydro.day(exp.hab.month$min.date,"Fed"))

exp.hab.month=exp.hab.month[order(exp.hab.month$min.date),]

plot(N.val~min.date,exp.hab.month,type="b")

n.samps.month=ddply(hab.dat2,c("CY","month"),summarise,min.date=min(date),N.val=N.obs(CELLCOUNT))

ylim.val=c(0,2000);by.y=500;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=date.fun(c("1953-10-01","2020-10-10"));xmaj=seq(xlim.val[1],xlim.val[2],"15 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years") #xmin=seq(xlim.val[1],xlim.val[2],"6 months")
# png(filename=paste0(plot.path,"HAB_sample_month.png"),width=7,height=2.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,3,0.5,1.5),oma=c(1,1,0.25,0.25));
plot(N.val~min.date,n.samps.month,type="n",ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,yaxs="i")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
polygon(c(expanded.dates,rev(expanded.dates)),
        c(-10,-10,ylim.val[2]+ylim.val[2]*0.5,ylim.val[2]+ylim.val[2]*0.5),col=adjustcolor("red",0.25),border=NA)
with(n.samps.month,segments(min.date,rep(0,length(min.date)),min.date,N.val))
axis_fun(1,xmaj,xmin,format(xmaj,"%m-%Y"),line=-0.5)
axis_fun(2,ymaj,ymaj,ymaj);box(lwd=1)
mtext(side=2,line=2.5,"Samples per Month")
mtext(side=1,line=1.75,"Date (Month-Year)")
dev.off()


# Test Train split
set.seed(111)
tr.index=sample(1:nrow(exp.hab.month),nrow(exp.hab.month)*0.7)
# tr.index=1:nrow(exp.hab.month);# No test train split
# tr.index=as.numeric(exp.hab.month$WY%in%seq(2013,2017,1))

# m1 <- gam(log10.mean~s(month,bs="cc",k=12)+s(DY,bs="tp")+ti(month,DY,bs=c("cc","tp")),
#           data=exp.hab.month[tr.index,], method = 'REML')

# m1 ----------------------------------------------------------------------
m1 <- gam(log10.mean~
            s(DoY,bs="cc",k=10)+
            s(CY,bs="tp",k=7)+
            ti(DoY,CY,bs=c("cc","tp"),k=c(9,5)),
          data=exp.hab.month[tr.index,], method = 'REML',knots=list(DoY=c(10,365.5)))


summary(m1)
nvar=3;layout(matrix(1:nvar,1,nvar))
plot(m1,residuals=T,pch=21)
dev.off()

nvar=4;layout(matrix(1:nvar,1,nvar))
gam.check(m1)
dev.off()

shapiro.test(m1$residuals)
acf(m1$residuals)

draw(m1)
pdat=with(exp.hab.month,expand.grid(DoY=seq(min(DoY,na.rm=T),max(DoY,na.rm=T),0.25),
                                    CY=seq(min(CY,na.rm=T),max(CY,na.rm=T),0.25)))

p.m1=predict(m1,newdata=pdat,type="terms",se=T)
pdat$fit.DoY=p.m1$fit[,1]
pdat$se.DoY=p.m1$se[,1]
pdat$fit.CY=p.m1$fit[,2]
pdat$se.CY=p.m1$se[,2]
pdat$fit.DoYCY=p.m1$fit[,3]

df.res <- df.residual(m1)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
# crit.t <- qnorm((1 - 0.95) / 2, lower.tail = FALSE)
pdat <- transform(pdat,
                   DoY.UCI = fit.DoY + (crit.t * se.DoY),
                   DoY.LCI = fit.DoY - (crit.t * se.DoY),
                   CY.UCI = fit.CY + (crit.t * se.CY),
                   CY.LCI = fit.CY - (crit.t * se.CY))

pred.org=predict(m1,type="terms")
partial.resids<-pred.org+residuals(m1)

hist(partial.resids[,1])
shapiro.test(partial.resids[,1])

hist(partial.resids[,2])
shapiro.test(partial.resids[,2])

hist(partial.resids[,3])
shapiro.test(partial.resids[,3])

# par(family="serif",mar=c(1,4,0.25,0.25),oma=c(1.75,2,0.75,0.5));
# layout(matrix(1:2,2,1))
# plot(m1,select=2,residuals=T,pch=21,ylim=c(-3,4))
# 
# pdat=pdat[order(pdat$CY,pdat$DoY),]
# plot(fit.CY~CY,pdat,ylim=c(-3,4),type="n")
# # http://zevross.com/blog/2014/09/15/recreate-the-gam-partial-regression-smooth-plots-from-r-package-mgcv-with-a-little-style/
# lines(fit.CY~CY,pdat)
# lines(CY.UCI~CY,pdat,lty=2)
# lines(CY.LCI~CY,pdat,lty=2)
# pred.org=predict(m1,type="terms")
# partial.resids<-pred.org+residuals(m1)
# points(m1$model$CY,partial.resids[,2],pch=21,bg="red")
# 
# par(family="serif",mar=c(1,4,0.25,0.25),oma=c(1.75,2,0.75,0.5));
# layout(matrix(1:2,2,1))
# plot(m1,select=1,residuals=T,pch=21,ylim=c(-3,4))
# 
# pdat=pdat[order(pdat$DoY,pdat$CY),]
# plot(fit.DoY~DoY,pdat,ylim=c(-3,4),type="n")
# lines(fit.DoY~DoY,pdat)
# lines(DoY.UCI~DoY,pdat,lty=2)
# lines(DoY.LCI~DoY,pdat,lty=2)
# pred.org=predict(m1,type="terms")
# partial.resids<-pred.org+residuals(m1)
# points(m1$model$DoY,partial.resids[,1],pch=21,bg="red")
# 
# par(family="serif",mar=c(1,4,0.25,0.25),oma=c(1.75,2,0.75,0.5));
# layout(matrix(1:2,2,1))
# plot(m1,select=3)
# 
# ylim.val=c(2012,2018);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2);ylen=length(ymaj)
# xlim.val=c(0,365);by.x=90;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2);xlen=length(xmaj)
# tmp.ma1=with(pdat,matrix(fit.DoYCY,ncol=length(unique(pdat$CY)),nrow=length(unique(pdat$DoY))))
# brk=50
# breaks.val=classIntervals(pdat$fit.DoYCY,style="equal",n=brk)
# pal=hcl.colors(n=brk,alpha=0.75)
# image(x=unique(pdat$DoY),y=unique(pdat$CY),z=tmp.ma1,
#       breaks=breaks.val$brks,col=pal,
#       ylim=ylim.val,xlim=xlim.val,axes=F,ann=F)
# contour(x=unique(pdat$DoY),y=unique(pdat$CY),z=tmp.ma1,add=T,drawlabels=F,lwd=1.25)
# axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
# mtext(side=3,adj=0,"ti(DoY,Year)")
# mtext(side=1,line=2.25,"DoY")
# mtext(side=2,line=3,"Year")


plot(exp.hab.month$log10.mean,type="b")

m.pred=exp.hab.month
m.pred=cbind(m.pred,data.frame(predict(m1,exp.hab.month,se.fit=T)))
m.pred$upr=with(m.pred, fit + (2*se.fit))
m.pred$lwr=with(m.pred, fit - (2*se.fit))
lines(m.pred$fit,col="Red")

ylim.val=c(0,7);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=expanded.dates;xmaj=seq(xlim.val[1],xlim.val[2],"1 years");xmin=seq(xlim.val[1],xlim.val[2],"6 months")
# png(filename=paste0(plot.path,"GAM_m1_HAB_month.png"),width=7,height=2.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,4,0.25,0.25),oma=c(1.75,2,0.75,0.5));

plot(log10.mean~min.date,exp.hab.month,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(exp.hab.month,lines(date.monCY,log10.mean,lwd=2))
with(exp.hab.month,points(date.monCY,log10.mean,pch=19))
with(m.pred,shaded.range(exp.hab.month$date.monCY,m.pred$upr,m.pred$lwr,bg="red",lty=0))
with(m.pred,lines(date.monCY,fit,col=adjustcolor("red",0.5),lwd=4))
axis_fun(1,xmaj,xmin,format(xmaj,"%m-%Y"),line=-0.5)
axis_fun(2,ymaj,ymaj,format(c(0,10^(ymaj[2:8])),scientific = F));box(lwd=1)
mtext(side=2,line=3.75,expression(paste(italic("K. brevis")," (cells L"^"-1",")")))
mtext(side=1,line=1.75,"Date (Month-Year)")
dev.off()

# png(filename=paste0(plot.path,"GAM_m1_draw_base.png"),width=8,height=2,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,3,1,1.5),oma=c(2,1,0.25,0.25));
layout(matrix(1:4,1,4),widths=c(1,1,1,0.5))

ylim.val=c(-4,4);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,365);by.x=90;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
pdat=pdat[order(pdat$DoY,pdat$CY),]
plot(fit.DoY~DoY,pdat,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
# lines(DoY.UCI~DoY,pdat,lty=2)
# lines(DoY.LCI~DoY,pdat,lty=2)
with(pdat,shaded.range(DoY,DoY.LCI,DoY.UCI,"grey",lty=1))
points(m1$model$DoY,partial.resids[,1],pch=19,col=adjustcolor("dodgerblue1",0.5))
lines(fit.DoY~DoY,pdat,lwd=2)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=3,adj=0,"s(DoY)")
mtext(side=1,line=2,"DoY")
mtext(side=2,line=2,"Effect")

ylim.val=c(-4,7);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(2012,2018);by.x=2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
pdat=pdat[order(pdat$CY,pdat$DoY),]
plot(fit.CY~CY,pdat,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(pdat,shaded.range(CY,CY.LCI,CY.UCI,"grey",lty=1))
points(m1$model$CY,partial.resids[,2],pch=19,col=adjustcolor("dodgerblue1",0.5))
lines(fit.CY~CY,pdat,lwd=2)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=3,adj=0,"s(Year)")
mtext(side=1,line=2,"Year")
mtext(side=2,line=2,"Effect")

ylim.val=c(2012,2018);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2);ylen=length(ymaj)
xlim.val=c(0,365);by.x=90;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2);xlen=length(xmaj)
tmp.ma1=with(pdat,matrix(fit.DoYCY,ncol=length(unique(pdat$CY)),nrow=length(unique(pdat$DoY))))
brk=50
breaks.val=classIntervals(pdat$fit.DoYCY,style="equal",n=brk)
pal=hcl.colors(n=brk,alpha=0.75)
image(x=unique(pdat$DoY),y=unique(pdat$CY),z=tmp.ma1,
      breaks=breaks.val$brks,col=pal,
      ylim=ylim.val,xlim=xlim.val,axes=F,ann=F)
contour(x=unique(pdat$DoY),y=unique(pdat$CY),z=tmp.ma1,add=T,drawlabels=F,lwd=1.25)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=3,adj=0,"ti(DoY,Year)")
mtext(side=1,line=2.25,"DoY")
mtext(side=2,line=3,"Year")

legend_image=as.raster(matrix(rev(pal),ncol=1))
par(xpd=NA,mar=c(2,1,1,0))
plot(c(0,1),c(0,1),type = 'n', axes = F,ann=F)
rasterImage(legend_image, 0, 0.25, 0.3,0.75)
text(x=0.3, y = seq(0.25,0.75,length.out=2), labels = format(round(range(breaks.val$brks),1)),cex=1,pos=4)
text(0.15,0.76,"Effect",pos=3,xpd=NA)
dev.off()

#
# theme_set(theme_minimal(base_size = 10, base_family = 'serif'))
# m1_draw=draw(m1,rug=F,nrow=1,residuals=T)
# ggsave(paste0(plot.path,"GAM_m1_draw.png"),m1_draw,device="png",height=2.5,width=9.25,unit="in")

cols="grey"#rainbow(length(mod.TP.all$fitted.values))

# tiff(filename=paste0(plot.path,"GAM_m1_diag.tiff"),width=6,height=5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
# png(filename=paste0(plot.path,"GAM_m1_diag.png"),width=6,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,3,1,0.5),oma=c(2,1,0.25,0.75));
layout(matrix(1:4,2,2))

ylim.val=c(-4,4);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(-4,4);by.x=2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
rstd=m1$residuals
qq.x=qq.function(m1$residuals)
plot(rstd~qq.x,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey80")
points(qq.x,rstd,pch=21,bg=adjustcolor(cols,0.5),col="grey",lwd=0.1,cex=1.25)
abline(0,1,lty=3)
axis_fun(1,line=-0.5,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2,"Deviance Residuals")
mtext(side=1,line=1.75,"Theoretical Quantiles")

ylim.val=c(-4,4);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,6);by.x=2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(m1$residuals~m1$linear.predictors,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey80")
abline(h=0)
with(m1,points(linear.predictors,residuals,pch=21,bg=adjustcolor(cols,0.5),col="grey",lwd=0.1,cex=1.25))
axis_fun(1,line=-0.5,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2,"Residuals")
mtext(side=1,line=1.75,"Linear Predictor")

hist.val=hist(m1$residuals,xlim=ylim.val,yaxs="i",col="grey",main=NA,plot=F)
xlim.val=ylim.val;by.x=by.y;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,50);by.y=10;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
hist(m1$residuals,xlim=xlim.val,ylim=ylim.val,yaxs="i",col="grey",main=NA,axes=F,ann=F)
axis_fun(1,xmaj,xmin,xmaj)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2,"Frequency")
mtext(side=1,line=1.75,"Residuals")

ylim.val=c(0,7);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=ylim.val;by.x=2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(m1$model$log10.mean~m1$fitted.values,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey80")
abline(0,1)
points(m1$fitted.values,m1$model$log10.mean,pch=21,bg=adjustcolor(cols,0.5),col="grey",lwd=0.1,cex=1.25)
axis_fun(1,xmaj,xmin,xmaj)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2,"Response")
mtext(side=1,line=1.75,"Fitted Values")
dev.off()

resid.val=m1$residuals
acf.m1.rslt=data.frame()
for(h in 0:24){
  #demean
  #x=sweep(as.matrix(wq.dat.xtab.mon$mean.TP),2,colMeans(as.matrix(wq.dat.xtab.mon$mean.TP),na.rm=T))
  lagged=lag(as.zoo(resid.val),-h,na.pad=T)
  tmp.dat=as.zoo(resid.val)
  stat=with(data.frame(lag=lagged,dat=tmp.dat),cor.test(lag,dat,method="pearson"))
  acf.m1.rslt=rbind(acf.m1.rslt,data.frame(lag=h,estimate=as.numeric(stat$estimate),pval=stat$p.value))
}

ylim.val=c(-1,1.1);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,24);by.x=5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
# png(filename=paste0(plot.path,"GAM_m1_Resid_ACF.png"),width=4,height=3.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.75,0.75),oma=c(2,1,0.5,0.5));
plot(estimate~lag,acf.m1.rslt,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
ci.val=qnorm((1+0.95)/2)/sqrt(length(m1$residuals))
#abline(h=c(ci.val,-ci.val),lty=2,lwd=1.5,col="blue")
polygon(c(-1,25,25,-1),c(ci.val,ci.val,-ci.val,-ci.val),col=adjustcolor("grey",0.5),border=0)
with(acf.m1.rslt,segments(lag,0,lag,estimate,lwd=1.5,lty=2))
with(acf.m1.rslt,points(lag,estimate,pch=21,bg=ifelse(pval<0.05,"indianred1","dodgerblue1"),lwd=0.01))
axis_fun(1,line=-0.5,xmaj,xmin,xmaj)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,expression(paste("ACF ",italic("r")["Pearson"])))
mtext(side=1,line=1.5,"Lag")
legend("topright",legend=c("\u03C1 < 0.05","\u03C1 > 0.05","95% CI"),
       pch=c(21,21,22),pt.bg=c("indianred1","dodgerblue1",adjustcolor("grey",0.5)),col=c("black","black",NA),
       lty=NA,lwd=c(0.1,0.1,0),pt.cex=1.5,cex=0.7,ncol=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()



# mse
m1.mse=mean((predict(m1,exp.hab.month[-tr.index,])-exp.hab.month[-tr.index,"log10.mean"])^2)
#rmse
m1.rmse=sqrt(mean((exp.hab.month[-tr.index,"log10.mean"]-predict(m1,exp.hab.month[-tr.index,]))^2))


# test
mod.pred=predict(m1,exp.hab.month[-tr.index,])
actuals_preds <-data.frame(cbind(actuals=exp.hab.month[-tr.index,"log10.mean"],predicted=mod.pred))
# Kling-Gupta Efficiency
r.val=with(actuals_preds,cor(predicted,actuals,method="pearson"))
alpha.val=with(actuals_preds,sd(predicted)/sd(actuals))
beta.val=with(actuals_preds,mean(predicted)/mean(actuals))
m1.KG=1-sqrt((r.val-1)^2 + (alpha.val-1)^2 + (beta.val-1)^2)


notidy_as_flextable_gam(m1)%>%font(fontname="TimesNewRoman",part="all")%>%print("pptx")

mod.sum.m1=notidy_tidy_gam(m1)
mod.sum.m1$model="m1"

mod.est.m1=notidy_glance_gam(m1)
mod.est.m1$model="m1"

## 
# m2 <- gam(log10.mean~s(month,bs="cc",k=12)+s(DY)+s(S79)+ti(month,DY,bs=c("cc","tp")),
#           data=exp.hab.month[tr.index,], method = 'REML')

# m2 ----------------------------------------------------------------------
m2 <- gam(log10.mean~s(DoY,bs="cc",k=12)+s(CY,bs="tp",k=7)+s(S79)+ti(DoY,CY,bs=c("cc","tp")),
          data=exp.hab.month[tr.index,], method = 'REML',knots=list(DoY=c(10,350)))
summary(m2)
nvar=4;layout(matrix(1:nvar,1,nvar))
plot(m2,residuals=T,pch=21)
dev.off()
nvar=4;layout(matrix(1:nvar,1,nvar))
gam.check(m2)

dev.off()

shapiro.test(m2$residuals)
acf(m2$residuals)

plot(exp.hab.month$log10.mean,type="b",ylim=c(0,10))
m.pred=exp.hab.month
m.pred=cbind(m.pred,data.frame(predict(m2,exp.hab.month,se.fit=T)))
m.pred$upr=with(m.pred, fit + (2*se.fit))
m.pred$lwr=with(m.pred, fit - (2*se.fit))
lines(m.pred$fit,col="Red")

ylim.val=c(0,7);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=expanded.dates;xmaj=seq(xlim.val[1],xlim.val[2],"1 years");xmin=seq(xlim.val[1],xlim.val[2],"6 months")
# png(filename=paste0(plot.path,"GAM_m2_HAB_month.png"),width=7,height=2.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,4,0.25,0.25),oma=c(1.75,2,0.75,0.5));

plot(log10.mean~min.date,exp.hab.month,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(exp.hab.month,lines(min.date,log10.mean,lwd=2))
with(exp.hab.month,points(min.date,log10.mean,pch=19))
with(m.pred,shaded.range(exp.hab.month$date.monCY,m.pred$upr,m.pred$lwr,bg="red",lty=0))
with(m.pred,lines(date.monCY,fit,col=adjustcolor("red",0.5),lwd=4))
axis_fun(1,xmaj,xmin,format(xmaj,"%m-%Y"),line=-0.5)
axis_fun(2,ymaj,ymaj,format(c(0,10^(ymaj[2:8])),scientific = F));box(lwd=1)
mtext(side=2,line=3.75,expression(paste(italic("K. brevis")," (cells L"^"-1",")")))
mtext(side=1,line=1.75,"Date (Month-Year)")
dev.off()


pdat2=with(exp.hab.month,expand.grid(DoY=seq(min(DoY,na.rm=T),max(DoY,na.rm=T),0.25),
                                     CY=seq(min(CY,na.rm=T),max(CY,na.rm=T),0.25),
                                     S79=seq(0,350000,5000)))

p.m2=predict(m2,newdata=pdat2,type="terms",se=T)
pdat2$fit.DoY=p.m2$fit[,1]
pdat2$se.DoY=p.m2$se[,1]
pdat2$fit.CY=p.m2$fit[,2]
pdat2$se.CY=p.m2$se[,2]
pdat2$fit.S79=p.m2$fit[,3]
pdat2$se.S79=p.m2$se[,3]
pdat2$fit.DoYCY=p.m2$fit[,4]

df.res <- df.residual(m2)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
# crit.t <- qnorm((1 - 0.95) / 2, lower.tail = FALSE)
pdat2 <- transform(pdat2,
                  DoY.UCI = fit.DoY + (crit.t * se.DoY),
                  DoY.LCI = fit.DoY - (crit.t * se.DoY),
                  CY.UCI = fit.CY + (crit.t * se.CY),
                  CY.LCI = fit.CY - (crit.t * se.CY),
                  S79.UCI = fit.S79 + (crit.t * se.S79),
                  S79.LCI = fit.S79 - (crit.t * se.S79))

pred.org=predict(m2,type="terms")
partial.resids<-pred.org+residuals(m2)

hist(partial.resids[,1])
shapiro.test(partial.resids[,1])

hist(partial.resids[,2])
shapiro.test(partial.resids[,2])

hist(partial.resids[,3])
shapiro.test(partial.resids[,3])

hist(partial.resids[,4])
shapiro.test(partial.resids[,4])


# png(filename=paste0(plot.path,"GAM_m2_draw_base.png"),width=8,height=2,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,3,1,1.5),oma=c(2,1,0.5,0.25));
layout(matrix(1:5,1,5),widths=c(1,1,1,1,0.5))

ylim.val=c(-4,4);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,365);by.x=90;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
pdat2=pdat2[order(pdat2$DoY,pdat2$CY),]
plot(fit.DoY~DoY,pdat2,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(pdat2,shaded.range(DoY,DoY.LCI,DoY.UCI,"grey",lty=1))
points(m2$model$DoY,partial.resids[,1],pch=19,col=adjustcolor("dodgerblue1",0.5))
lines(fit.DoY~DoY,pdat2,lwd=2)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=3,adj=0,"s(DoY)")
mtext(side=1,line=2,"DoY")
mtext(side=2,line=2,"Effect")

ylim.val=c(-4,7);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(2012,2018);by.x=2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
pdat2=pdat2[order(pdat2$CY,pdat2$DoY),]
plot(fit.CY~CY,pdat2,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(pdat2,shaded.range(CY,CY.LCI,CY.UCI,"grey",lty=1))
points(m2$model$CY,partial.resids[,2],pch=19,col=adjustcolor("dodgerblue1",0.5))
lines(fit.CY~CY,pdat2,lwd=2)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=3,adj=0,"s(Year)")
mtext(side=1,line=2,"Year")
mtext(side=2,line=2,"Effect")

ylim.val=c(-4,4);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,350)*10e2;by.x=100*10e2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
pdat2=pdat2[order(pdat2$S79),]
plot(fit.S79~S79,pdat2,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(pdat2,shaded.range(S79,S79.LCI,S79.UCI,"grey",lty=1))
points(m2$model$S79,partial.resids[,3],pch=19,col=adjustcolor("dodgerblue1",0.5))
lines(fit.S79~S79,pdat2,lwd=2)
axis_fun(1,xmaj,xmin,xmaj/10e2,line=-0.5);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=3,adj=0,expression(paste("s(Q"["S79"],")")))
mtext(side=1,line=2,expression(paste("Q"["S79"]," (x1000 cfs)")))
mtext(side=2,line=2,"Effect")

ylim.val=c(2012,2018);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2);ylen=length(ymaj)
xlim.val=c(0,365);by.x=90;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2);xlen=length(xmaj)
tmp.ma1=with(pdat2,matrix(fit.DoYCY,ncol=length(unique(pdat2$CY)),nrow=length(unique(pdat2$DoY))))
brk=50
breaks.val=classIntervals(pdat2$fit.DoYCY,style="equal",n=brk)
pal=hcl.colors(n=brk,alpha=0.75)
image(x=unique(pdat2$DoY),y=unique(pdat2$CY),z=tmp.ma1,
      breaks=breaks.val$brks,col=pal,
      ylim=ylim.val,xlim=xlim.val,axes=F,ann=F)
contour(x=unique(pdat2$DoY),y=unique(pdat2$CY),z=tmp.ma1,add=T,drawlabels=F,lwd=1.25)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=3,adj=0,"ti(DoY,Year)")
mtext(side=1,line=2.25,"DoY")
mtext(side=2,line=3,"Year")

legend_image=as.raster(matrix(rev(pal),ncol=1))
par(xpd=NA,mar=c(2,1,1,0))
plot(c(0,1),c(0,1),type = 'n', axes = F,ann=F)
rasterImage(legend_image, 0, 0.25, 0.3,0.75)
text(x=0.3, y = seq(0.25,0.75,length.out=2), labels = format(round(range(breaks.val$brks),1)),cex=1,pos=4)
text(0.15,0.76,"Effect",pos=3,xpd=NA)
dev.off()



# m2_draw=draw(m2,rug=F,nrow=1,residuals=T)
# ggsave(paste0(plot.path,"GAM_m2_draw.png"),m2_draw,device="png",height=2.5,width=9.25,unit="in")


cols="grey"#rainbow(length(mod.TP.all$fitted.values))

# tiff(filename=paste0(plot.path,"GAM_m2_diag.tiff"),width=6,height=5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
# png(filename=paste0(plot.path,"GAM_m2_diag.png"),width=6,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,3,1,0.5),oma=c(2,1,0.25,0.75));
layout(matrix(1:4,2,2))

ylim.val=c(-3,3);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(-3,3);by.x=1;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
rstd=m2$residuals
qq.x=qq.function(m2$residuals)
plot(rstd~qq.x,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey80")
points(qq.x,rstd,pch=21,bg=adjustcolor(cols,0.5),col="grey",lwd=0.1,cex=1.25)
abline(0,1,lty=3)
axis_fun(1,line=-0.5,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2,"Deviance Residuals")
mtext(side=1,line=1.75,"Theoretical Quantiles")

ylim.val=c(-3,3);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,6);by.x=2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(m2$residuals~m2$linear.predictors,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey80")
abline(h=0)
with(m2,points(linear.predictors,residuals,pch=21,bg=adjustcolor(cols,0.5),col="grey",lwd=0.1,cex=1.25))
axis_fun(1,line=-0.5,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2,"Residuals")
mtext(side=1,line=1.75,"Linear Predictor")

hist.val=hist(m2$residuals,plot=F)
xlim.val=ylim.val;by.x=by.y;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,20);by.y=5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
hist(m2$residuals,xlim=xlim.val,ylim=ylim.val,yaxs="i",col="grey",main=NA,axes=F,ann=F)
axis_fun(1,xmaj,xmin,xmaj)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2,"Frequency")
mtext(side=1,line=1.75,"Residuals")

ylim.val=c(0,7);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=ylim.val;by.x=2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(m2$model$log10.mean~m2$fitted.values,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey80")
abline(0,1)
points(m2$fitted.values,m2$model$log10.mean,pch=21,bg=adjustcolor(cols,0.5),col="grey",lwd=0.1,cex=1.25)
axis_fun(1,xmaj,xmin,xmaj)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2,"Response")
mtext(side=1,line=1.75,"Fitted Values")
dev.off()

resid.val=m2$residuals
acf.m2.rslt=data.frame()
for(h in 0:24){
  #demean
  #x=sweep(as.matrix(wq.dat.xtab.mon$mean.TP),2,colMeans(as.matrix(wq.dat.xtab.mon$mean.TP),na.rm=T))
  lagged=lag(as.zoo(resid.val),-h,na.pad=T)
  tmp.dat=as.zoo(resid.val)
  stat=with(data.frame(lag=lagged,dat=tmp.dat),cor.test(lag,dat,method="pearson"))
  acf.m2.rslt=rbind(acf.m2.rslt,data.frame(lag=h,estimate=as.numeric(stat$estimate),pval=stat$p.value))
}

ylim.val=c(-1,1.1);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,24);by.x=5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
# png(filename=paste0(plot.path,"GAM_m2_Resid_ACF.png"),width=4,height=3.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.75,0.75),oma=c(2,1,0.5,0.5));
plot(estimate~lag,acf.m2.rslt,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
ci.val=qnorm((1+0.95)/2)/sqrt(length(m2$residuals))
#abline(h=c(ci.val,-ci.val),lty=2,lwd=1.5,col="blue")
polygon(c(-1,25,25,-1),c(ci.val,ci.val,-ci.val,-ci.val),col=adjustcolor("grey",0.5),border=0)
with(acf.m2.rslt,segments(lag,0,lag,estimate,lwd=1.5,lty=2))
with(acf.m2.rslt,points(lag,estimate,pch=21,bg=ifelse(pval<0.05,"indianred1","dodgerblue1"),lwd=0.01))
axis_fun(1,line=-0.5,xmaj,xmin,xmaj)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,expression(paste("ACF ",italic("r")["Pearson"])))
mtext(side=1,line=1.5,"Lag")
legend("topright",legend=c("\u03C1 < 0.05","\u03C1 > 0.05","95% CI"),
       pch=c(21,21,22),pt.bg=c("indianred1","dodgerblue1",adjustcolor("grey",0.5)),col=c("black","black",NA),
       lty=NA,lwd=c(0.1,0.1,0),pt.cex=1.5,cex=0.7,ncol=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()

#mse
m2.mse=mean((exp.hab.month[-tr.index,"log10.mean"]-predict(m2,exp.hab.month[-tr.index,]))^2)
#rmse
m2.rmse=sqrt(mean((exp.hab.month[-tr.index,"log10.mean"]-predict(m2,exp.hab.month[-tr.index,]))^2))

# test
mod.pred=predict(m2,exp.hab.month[-tr.index,])
actuals_preds <-data.frame(cbind(actuals=exp.hab.month[-tr.index,"log10.mean"],predicted=mod.pred))
# Kling-Gupta Efficiency
r.val=with(actuals_preds,cor(predicted,actuals,method="pearson"))
alpha.val=with(actuals_preds,sd(predicted)/sd(actuals))
beta.val=with(actuals_preds,mean(predicted)/mean(actuals))
m2.KG=1-sqrt((r.val-1)^2 + (alpha.val-1)^2 + (beta.val-1)^2)

notidy_as_flextable_gam(m2)%>%font(fontname="TimesNewRoman",part="all")%>%print("pptx")

mod.sum.m2=notidy_tidy_gam(m2)
mod.sum.m2$model="m2"

mod.est.m2=notidy_glance_gam(m2)
mod.est.m2$model="m2"


##
## period of change
## see https://fromthebottomoftheheap.net/2014/05/15/identifying-periods-of-change-with-gams/
tmpf <- tempfile()
download.file("https://gist.github.com/gavinsimpson/e73f011fdaaab4bb5a30/raw/82118ee30c9ef1254795d2ec6d356a664cc138ab/Deriv.R",
              tmpf)
source(tmpf)

exp.hab.month$DoY
exp.hab.month$CY
exp.hab.month$S79

pdat2=with(exp.hab.month,expand.grid(DoY=seq(min(DoY,na.rm=T),max(DoY,na.rm=T),0.25),
                                       CY=seq(min(CY,na.rm=T),max(CY,na.rm=T),0.25),
                                       S79=seq(0,350000,5000)))

p2 <- predict(m2, newdata=pdat2,type = "terms", se.fit = TRUE)
pdat2$p2_S79=p2$fit[,3]
pdat2$se2_S79 = p2$se.fit[,3]

df.res <- df.residual(m2)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
pdat2 <- transform(pdat2,
                     upper = p2_S79 + (crit.t * se2_S79),
                     lower = p2_S79 - (crit.t * se2_S79))
pdat2=pdat2[order(pdat2$CY,pdat2$DoY),]

m2.d <- derivatives(m2,newdata=pdat2,term='s(S79)',type = "central",interval="confidence",ncores=2)

m2.dsig <- signifD(pdat2$p2,
                     d=m2.d$derivative,
                     m2.d$upper,m2.d$lower)
pdat2$dsig.incr=unlist(m2.dsig$incr)
pdat2$dsig.decr=unlist(m2.dsig$decr)

subset(pdat2,is.na(dsig.incr)==F)$S79
subset(pdat2,is.na(dsig.decr)==F)$S79


###

# m3 <- gam(log10.mean~s(month,bs="cc",k=12)+s(DY)+s(mean.stg)+ti(month,DY,bs=c("cc","tp")),
#           data=exp.hab.month[tr.index,], method = 'REML')

# m3 ----------------------------------------------------------------------
m3 <- gam(log10.mean~s(DoY,bs="cc",k=12)+s(CY,bs="tp",k=7)+s(mean.stg)+ti(DoY,CY,bs=c("cc","tp")),
          data=exp.hab.month[tr.index,], method = 'REML',knots=list(DoY=c(10,350)))

summary(m3)
nvar=4;layout(matrix(1:nvar,1,nvar))
plot(m3,residuals=T,pch=21)
dev.off()
nvar=4;layout(matrix(1:nvar,1,nvar))
gam.check(m3)

dev.off()

shapiro.test(m3$residuals)
acf(m3$residuals)

plot(exp.hab.month$log10.mean,type="b",ylim=c(0,10))
m.pred=exp.hab.month
m.pred=cbind(m.pred,data.frame(predict(m3,exp.hab.month,se.fit=T)))
m.pred$upr=with(m.pred, fit + (2*se.fit))
m.pred$lwr=with(m.pred, fit - (2*se.fit))
lines(m.pred$fit,col="Red")

ylim.val=c(0,7);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=expanded.dates;xmaj=seq(xlim.val[1],xlim.val[2],"1 years");xmin=seq(xlim.val[1],xlim.val[2],"6 months")
# png(filename=paste0(plot.path,"GAM_m3_HAB_month.png"),width=7,height=2.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,4,0.25,0.25),oma=c(1.75,2,0.75,0.5));

plot(log10.mean~min.date,exp.hab.month,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(exp.hab.month,lines(min.date,log10.mean,lwd=2))
with(exp.hab.month,points(min.date,log10.mean,pch=19))
with(m.pred,shaded.range(exp.hab.month$date.monCY,m.pred$upr,m.pred$lwr,bg="red",lty=0))
with(m.pred,lines(date.monCY,fit,col=adjustcolor("red",0.5),lwd=4))
axis_fun(1,xmaj,xmin,format(xmaj,"%m-%Y"),line=-0.5)
axis_fun(2,ymaj,ymaj,format(c(0,10^(ymaj[2:8])),scientific = F));box(lwd=1)
mtext(side=2,line=3.75,expression(paste(italic("K. brevis")," (cells L"^"-1",")")))
mtext(side=1,line=1.75,"Date (Month-Year)")
dev.off()

pdat3=with(exp.hab.month,expand.grid(DoY=seq(min(DoY,na.rm=T),max(DoY,na.rm=T),0.25),
                                     CY=seq(min(CY,na.rm=T),max(CY,na.rm=T),0.25),
                                     mean.stg=seq(min(mean.stg,na.rm=T),max(mean.stg,na.rm=T),0.1)))

p.m3=predict(m3,newdata=pdat3,type="terms",se=T)
pdat3$fit.DoY=p.m3$fit[,1]
pdat3$se.DoY=p.m3$se[,1]
pdat3$fit.CY=p.m3$fit[,2]
pdat3$se.CY=p.m3$se[,2]
pdat3$fit.stage=p.m3$fit[,3]
pdat3$se.stage=p.m3$se[,3]
pdat3$fit.DoYCY=p.m3$fit[,4]

df.res <- df.residual(m3)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
# crit.t <- qnorm((1 - 0.95) / 2, lower.tail = FALSE)
pdat3 <- transform(pdat3,
                   DoY.UCI = fit.DoY + (crit.t * se.DoY),
                   DoY.LCI = fit.DoY - (crit.t * se.DoY),
                   CY.UCI = fit.CY + (crit.t * se.CY),
                   CY.LCI = fit.CY - (crit.t * se.CY),
                   stage.UCI = fit.stage + (crit.t * se.stage),
                   stage.LCI = fit.stage - (crit.t * se.stage))

pred.org=predict(m3,type="terms")
partial.resids<-pred.org+residuals(m3)

hist(partial.resids[,1])
shapiro.test(partial.resids[,1])

hist(partial.resids[,2])
shapiro.test(partial.resids[,2])

hist(partial.resids[,3])
shapiro.test(partial.resids[,3])

hist(partial.resids[,4])
shapiro.test(partial.resids[,4])


# png(filename=paste0(plot.path,"GAM_m3_draw_base.png"),width=8,height=2,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,3,1,1.5),oma=c(2,1,0.5,0.25));
layout(matrix(1:5,1,5),widths=c(1,1,1,1,0.5))

ylim.val=c(-4,4);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,365);by.x=90;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
pdat3=pdat3[order(pdat3$DoY,pdat3$CY),]
plot(fit.DoY~DoY,pdat3,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(pdat3,shaded.range(DoY,DoY.LCI,DoY.UCI,"grey",lty=1))
points(m3$model$DoY,partial.resids[,1],pch=19,col=adjustcolor("dodgerblue1",0.5))
lines(fit.DoY~DoY,pdat3,lwd=2)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=3,adj=0,"s(DoY)")
mtext(side=1,line=2,"DoY")
mtext(side=2,line=2,"Effect")

ylim.val=c(-4,7);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(2012,2018);by.x=2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
pdat3=pdat3[order(pdat3$CY,pdat3$DoY),]
plot(fit.CY~CY,pdat3,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(pdat3,shaded.range(CY,CY.LCI,CY.UCI,"grey",lty=1))
points(m3$model$CY,partial.resids[,2],pch=19,col=adjustcolor("dodgerblue1",0.5))
lines(fit.CY~CY,pdat3,lwd=2)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=3,adj=0,"s(Year)")
mtext(side=1,line=2,"Year")
mtext(side=2,line=2,"Effect")

ylim.val=c(-4,4);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(11,17);by.x=2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
pdat3=pdat3[order(pdat3$mean.stg),]
plot(fit.stage~mean.stg,pdat3,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(pdat3,shaded.range(mean.stg,stage.LCI,stage.UCI,"grey",lty=1))
points(m3$model$mean.stg,partial.resids[,3],pch=19,col=adjustcolor("dodgerblue1",0.5))
lines(fit.stage~mean.stg,pdat3,lwd=2)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=3,adj=0,"s(mean Lake Stage)",cex=0.85)
mtext(side=1,line=2,"Stage (Ft, NGVD29)")
mtext(side=2,line=2,"Effect")

ylim.val=c(2012,2018);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2);ylen=length(ymaj)
xlim.val=c(0,365);by.x=90;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2);xlen=length(xmaj)
tmp.ma1=with(pdat3,matrix(fit.DoYCY,ncol=length(unique(pdat3$CY)),nrow=length(unique(pdat3$DoY))))
brk=50
breaks.val=classIntervals(pdat3$fit.DoYCY,style="equal",n=brk)
pal=hcl.colors(n=brk,alpha=0.75)
image(x=unique(pdat3$DoY),y=unique(pdat3$CY),z=tmp.ma1,
      breaks=breaks.val$brks,col=pal,
      ylim=ylim.val,xlim=xlim.val,axes=F,ann=F)
contour(x=unique(pdat3$DoY),y=unique(pdat3$CY),z=tmp.ma1,add=T,drawlabels=F,lwd=1.25)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=3,adj=0,"ti(DoY,Year)")
mtext(side=1,line=2.25,"DoY")
mtext(side=2,line=3,"Year")

legend_image=as.raster(matrix(rev(pal),ncol=1))
par(xpd=NA,mar=c(2,1,1,0))
plot(c(0,1),c(0,1),type = 'n', axes = F,ann=F)
rasterImage(legend_image, 0, 0.25, 0.3,0.75)
text(x=0.3, y = seq(0.25,0.75,length.out=2), labels = format(round(range(breaks.val$brks),1)),cex=1,pos=4)
text(0.15,0.76,"Effect",pos=3,xpd=NA)
dev.off()

# m3_draw=draw(m3,rug=F,nrow=1,residuals=T)
# ggsave(paste0(plot.path,"GAM_m3_draw.png"),m3_draw,device="png",height=2.5,width=9.25,unit="in")

cols="grey"#rainbow(length(mod.TP.all$fitted.values))

# tiff(filename=paste0(plot.path,"GAM_m3_diag.tiff"),width=6,height=5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
# png(filename=paste0(plot.path,"GAM_m3_diag.png"),width=6,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,3,1,0.5),oma=c(2,1,0.25,0.75));
layout(matrix(1:4,2,2))

ylim.val=c(-3,3);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(-3,3);by.x=1;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
rstd=m3$residuals
qq.x=qq.function(m3$residuals)
plot(rstd~qq.x,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey80")
points(qq.x,rstd,pch=21,bg=adjustcolor(cols,0.5),col="grey",lwd=0.1,cex=1.25)
abline(0,1,lty=3)
axis_fun(1,line=-0.5,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2,"Deviance Residuals")
mtext(side=1,line=1.75,"Theoretical Quantiles")

ylim.val=c(-3,3);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,6);by.x=2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(m3$residuals~m3$linear.predictors,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey80")
abline(h=0)
with(m3,points(linear.predictors,residuals,pch=21,bg=adjustcolor(cols,0.5),col="grey",lwd=0.1,cex=1.25))
axis_fun(1,line=-0.5,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2,"Residuals")
mtext(side=1,line=1.75,"Linear Predictor")

hist.val=hist(m3$residuals,plot=F)
xlim.val=ylim.val;by.x=by.y;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,15);by.y=3;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
hist(m3$residuals,xlim=xlim.val,ylim=ylim.val,yaxs="i",col="grey",main=NA,axes=F,ann=F)
axis_fun(1,xmaj,xmin,xmaj)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2,"Frequency")
mtext(side=1,line=1.75,"Residuals")

ylim.val=c(0,7);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=ylim.val;by.x=2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(m3$model$log10.mean~m3$fitted.values,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey80")
abline(0,1)
points(m3$fitted.values,m3$model$log10.mean,pch=21,bg=adjustcolor(cols,0.5),col="grey",lwd=0.1,cex=1.25)
axis_fun(1,xmaj,xmin,xmaj)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2,"Response")
mtext(side=1,line=1.75,"Fitted Values")
dev.off()

resid.val=m3$residuals
acf.m3.rslt=data.frame()
for(h in 0:24){
  #demean
  #x=sweep(as.matrix(wq.dat.xtab.mon$mean.TP),2,colMeans(as.matrix(wq.dat.xtab.mon$mean.TP),na.rm=T))
  lagged=lag(as.zoo(resid.val),-h,na.pad=T)
  tmp.dat=as.zoo(resid.val)
  stat=with(data.frame(lag=lagged,dat=tmp.dat),cor.test(lag,dat,method="pearson"))
  acf.m3.rslt=rbind(acf.m3.rslt,data.frame(lag=h,estimate=as.numeric(stat$estimate),pval=stat$p.value))
}

ylim.val=c(-1,1.1);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,24);by.x=5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
# png(filename=paste0(plot.path,"GAM_m3_Resid_ACF.png"),width=4,height=3.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.75,0.75),oma=c(2,1,0.5,0.5));
plot(estimate~lag,acf.m3.rslt,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
ci.val=qnorm((1+0.95)/2)/sqrt(length(m3$residuals))
#abline(h=c(ci.val,-ci.val),lty=2,lwd=1.5,col="blue")
polygon(c(-1,25,25,-1),c(ci.val,ci.val,-ci.val,-ci.val),col=adjustcolor("grey",0.5),border=0)
with(acf.m3.rslt,segments(lag,0,lag,estimate,lwd=1.5,lty=2))
with(acf.m3.rslt,points(lag,estimate,pch=21,bg=ifelse(pval<0.05,"indianred1","dodgerblue1"),lwd=0.01))
axis_fun(1,line=-0.5,xmaj,xmin,xmaj)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,expression(paste("ACF ",italic("r")["Pearson"])))
mtext(side=1,line=1.5,"Lag")
legend("topright",legend=c("\u03C1 < 0.05","\u03C1 > 0.05","95% CI"),
       pch=c(21,21,22),pt.bg=c("indianred1","dodgerblue1",adjustcolor("grey",0.5)),col=c("black","black",NA),
       lty=NA,lwd=c(0.1,0.1,0),pt.cex=1.5,cex=0.7,ncol=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()

# #mse
m3.mse=mean((exp.hab.month[-tr.index,"log10.mean"]-predict(m3,exp.hab.month[-tr.index,]))^2)
#rmse
m3.rmse=sqrt(mean((exp.hab.month[-tr.index,"log10.mean"]-predict(m3,exp.hab.month[-tr.index,]))^2))

# test
mod.pred=predict(m3,exp.hab.month[-tr.index,])
actuals_preds <-data.frame(cbind(actuals=exp.hab.month[-tr.index,"log10.mean"],predicted=mod.pred))
# Kling-Gupta Efficiency
r.val=with(actuals_preds,cor(predicted,actuals,method="pearson"))
alpha.val=with(actuals_preds,sd(predicted)/sd(actuals))
beta.val=with(actuals_preds,mean(predicted)/mean(actuals))
m3.KG=1-sqrt((r.val-1)^2 + (alpha.val-1)^2 + (beta.val-1)^2)


notidy_as_flextable_gam(m3)#%>%font(fontname="TimesNewRoman",part="all")%>%print("pptx")

mod.sum.m3=notidy_tidy_gam(m3)
mod.sum.m3$model="m3"

mod.est.m3=notidy_glance_gam(m3)
mod.est.m3$model="m3"

mod.sum.all=rbind(mod.sum.m1,mod.sum.m2,mod.sum.m3)
mod.est.all=rbind(mod.est.m1,mod.est.m2,mod.est.m3)
# write.csv(mod.sum.all,paste0(export.path,"gam_mod_sum.csv"),row.names = F)
# write.csv(mod.est.all,paste0(export.path,"gam_est_sum.csv"),row.names = F)

data.frame(model=c("Season","Season + Discharge","Season + Lake Stage"),
           adjR2=mod.est.all$adj.r.squared,
           dev.exp=mod.est.all$deviance,
           MSE=c(m1.mse,m2.mse,m3.mse),
           KG=c(m1.KG,m2.KG,m3.KG),
           df=mod.est.all$df,
           AIC=mod.est.all$AIC)%>%
  flextable()%>%
  colformat_double(j=2:5,digits=2)%>%
  colformat_double(j=6:7,digits=1)%>%
  bold(part="header")%>%
  set_header_labels("model"="Candidate Model",
                    "adjR2"="Adj R\u00B2",
                    "dev.exp"="Deviance\nExplained",
                    "MSE"="Mean Square\nError",
                    "df"="Degrees of\nFreedom")%>%
  align(j=2:6,align="center",part="all")%>%
  width(width=c(1.75,1,1,1.25,0.75,1,0.75))%>%
  footnote(i=1,j=3,part="header",
           value=as_paragraph(" Deviance Explained is similar to unadjusted R\u00B2"),
           ref_symbols = c(" 1"))%>%
  footnote(i=1,j=5,part="header",
           value=as_paragraph(" KG = Kling-Gupta coefficient"),
           ref_symbols = c(" 2"))%>%
  footnote(i=1,j=7,part="header",
           value=as_paragraph(" AIC = Akaike Information Criterion"),
           ref_symbols = c(" 3"))%>%
  font(fontname="TimesNewRoman",part="all")# %>%print("pptx")


## Period of Change
pdat3=with(exp.hab.month,expand.grid(DoY=seq(min(DoY,na.rm=T),max(DoY,na.rm=T),0.25),
                                    CY=seq(min(CY,na.rm=T),max(CY,na.rm=T),0.25),
                                    mean.stg=seq(min(mean.stg,na.rm=T),max(mean.stg,na.rm=T),0.1)))

p2 <- predict(m3, newdata=pdat3,type = "terms", se.fit = TRUE)
head(p2$fit)
pdat3$p2=p2$fit[,3]
pdat3$se2 = p2$se.fit[,3]

df.res <- df.residual(m3)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
pdat3 <- transform(pdat3,
                  upper = p2 + (crit.t * se2),
                  lower = p2 - (crit.t * se2))
pdat3=pdat3[order(pdat3$CY,pdat3$DoY),]

m3.d <- derivatives(m3,newdata=pdat3,term='s(mean.stg)',type = "central",interval="confidence",ncores=2)
m3.dsig <- signifD(pdat3$p2,
                   d=m3.d$derivative,
                   m3.d$upper,m3.d$lower)
pdat3$dsig.incr=unlist(m3.dsig$incr)
pdat3$dsig.decr=unlist(m3.dsig$decr)

subset(pdat3,is.na(dsig.incr)==F)
subset(pdat3,is.na(dsig.decr)==F)

range(subset(pdat3,is.na(dsig.incr)==F)$mean.stg)

pdat3=pdat3[order(pdat3$mean.stg),]

# png(filename=paste0(plot.path,"m3_stage_change.png"),width=5,height=3.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,2,0.5,0.5),oma=c(1,2,0.75,0.5));
ylim.val=c(-4,4);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(11,17);by.x=2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)

plot(p2 ~ mean.stg, data = pdat3,type="n",ylim=ylim.val,xlim=xlim.val,ann=F,axes=F)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
lines(p2 ~ mean.stg, data = pdat3)
lines(upper ~ mean.stg, data = pdat3, lty = "dashed")
lines(lower ~ mean.stg, data = pdat3, lty = "dashed")
lines(dsig.incr ~ mean.stg, data = pdat3, col = "red", lwd = 3,lty=1)
lines(dsig.decr ~ mean.stg, data = pdat3, col = "blue", lwd = 3,lty=1)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=1,line=1.5,"Stage (Ft, NGVD29)")
mtext(side=2,line=2.5,"Effect")
legend("topleft",legend=c("Fitted (mean Lake Stage)","95% CI (WY)","Significant Change"),
       lty=c(1,2,1,0),col=c("black","black","red"),
       pch=c(NA,NA,NA),pt.bg=c(NA,NA,NA),
       pt.cex=1,ncol=1,cex=0.6,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1)
dev.off()

exp.hab.month$Stg_chg=with(exp.hab.month,ifelse(mean.stg>min(subset(pdat3,is.na(dsig.incr)==F)$mean.stg)&
                                                  mean.stg<max(subset(pdat3,is.na(dsig.incr)==F)$mean.stg),1,0))

ylim.val=c(0.001,1000)*10e2
boxplot(mean.val~Stg_chg,exp.hab.month,outline=F,log="y",
        ylim=ylim.val,col=adjustcolor(c("indianred1","dodgerblue1"),0.5))
kruskal.test(mean.val~Stg_chg,exp.hab.month)
kruskal.test(log10.mean~Stg_chg,exp.hab.month)
## 


AIC(m1,m2,m3)
anova(m1,m2,test="Chisq")
anova(m1,m3,test="Chisq")
anova(m2,m3,test="Chisq")

mod.comp=anova(m1,m2,m3,test="F")
mod.comp=cbind(data.frame(model=c("Season","Season + Discharge","Season + Lake Stage")),mod.comp)

mod.comp%>%
  flextable()%>%
  bold(part="header")%>%
  colformat_double(j=2:7,digits=2,na_str=" ")%>%
  set_header_labels("model"="Candidate Model",
                    "Resid.Df"="Residual\nDOF",
                    "Resid.Dev"="Residual\nDeviance",
                    "Df"="DOF",
                    "Deviance"="Deviance",
                    "Fstat"="F-Statistic",
                    "Pr(>F)"="\u03C1-value")%>%
  align(j=2:7,align="center",part="all")%>%
  width(width=c(1.75,1,1,0.75,0.75,0.5,0.75))%>%
  font(fontname="TimesNewRoman",part="all")%>%print("pptx")
  

# test
mod.pred=predict(m2,exp.hab.month[-tr.index,])
actuals_preds <-data.frame(cbind(actuals=exp.hab.month[-tr.index,"log10.mean"],predicted=mod.pred))
plot(actuals~predicted,actuals_preds);abline(0,1)
cor.test(actuals_preds$actuals,actuals_preds$predicted,method="spearman")

# Nash-Sutcliffe model efficiency coefficient
with(actuals_preds,1-(sum((actuals-predicted)^2)/sum((actuals-mean(actuals))^2)))
# Kling-Gupta Efficiency
r.val=with(actuals_preds,cor(predicted,actuals,method="pearson"))
alpha.val=with(actuals_preds,sd(predicted)/sd(actuals))
beta.val=with(actuals_preds,mean(predicted)/mean(actuals))
1-sqrt((r.val-1)^2 + (alpha.val-1)^2 + (beta.val-1)^2)

visreg::visreg2d(m2, xvar='DoY', yvar='CY', scale='response')
visreg::visreg2d(m3, xvar='DoY', yvar='CY', scale='response')


#####

m4 <- gam(log10.mean~s(month,bs="cc",k=12)+s(DY)+s(S79_TN)+ti(month,DY,bs=c("cc","tp")),
          data=exp.hab.month[tr.index,], method = 'REML')
summary(m4)
nvar=5;layout(matrix(1:nvar,1,nvar))
plot(m4,residuals=T,pch=21)
dev.off()
nvar=4;layout(matrix(1:nvar,1,nvar))
gam.check(m4)

dev.off()

shapiro.test(m4$residuals)
acf(m4$residuals)

plot(exp.hab.month$log10.mean,type="b",ylim=c(0,10))
m.pred=data.frame(predict(m4,exp.hab.month,se.fit=T))
lines(m.pred$fit,col=adjustcolor("red",0.5),lwd=2)
lines(predict(m1,exp.hab.month),col="olivedrab3",lwd=2)
lines(predict(m2,exp.hab.month),col="dodgerblue1",lwd=2)
lines(predict(m3,exp.hab.month),col="orange",lwd=2)
#mse
# mean((exp.hab.month[-tr.index,"log10.mean"]-predict(m4,exp.hab.month[-tr.index,]))^2)

AIC(m1,m2,m3,m4)
AIC(m1,m2,m3)



##### 
ylim.val=c(0,7);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=expanded.dates;xmaj=seq(xlim.val[1],xlim.val[2],"1 years");xmin=seq(xlim.val[1],xlim.val[2],"6 months")
# png(filename=paste0(plot.path,"GAM_all_HAB_month.png"),width=7,height=3,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,4,0.25,0.25),oma=c(1.75,2,0.75,0.5));
layout(matrix(1:2,2,1),heights=c(1,0.5))
cols=c("indianred1","dodgerblue1","darkseagreen3")#hcl.colors(3)
plot(log10.mean~min.date,exp.hab.month,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(exp.hab.month,lines(min.date,log10.mean,lwd=2,col="grey"))
with(exp.hab.month,points(min.date,log10.mean,pch=19,col="grey"))
lines(exp.hab.month$min.date,predict(m1,exp.hab.month),col=cols[1],lwd=2)
lines(exp.hab.month$min.date,predict(m2,exp.hab.month),col=cols[2],lwd=2)
lines(exp.hab.month$min.date,predict(m3,exp.hab.month),col=cols[3],lwd=2)
axis_fun(1,xmaj,xmin,format(xmaj,"%m-%Y"),line=-0.5)
axis_fun(2,ymaj,ymaj,format(c(0,10^(ymaj[2:8])),scientific = F));box(lwd=1)
mtext(side=2,line=4,expression(paste(italic("K. brevis")," (cells L"^"-1",")")))
mtext(side=1,line=1.75,"Date (Month-Year)")

plot(0:1,0:1,type = 'n', axes = F,ann=F)
leg.text=c("Observed Data","Season","Season + Discharge","Season + Lake Stage")
legend(0.5,-0.25,legend=leg.text,
       pch=c(19,NA,NA,NA),
       lty=c(1,1,1,1),lwd=c(1,2,2,2),
       col=c("grey",cols),
       pt.cex=1,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)

dev.off()

