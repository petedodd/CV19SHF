## process sitrep
## new version using revised PHE case version
## should be run before modelling script
## p.j.dodd@sheffield.ac.uk
library(here)
library(readxl)
library(data.table)
library(lubridate)
library(ggpubr)
library(glue)
library(scales)
library(ggplot2)

getdate <- function()glue(gsub('-','_',Sys.Date()))
td <- getdate()                         #date stamp

## =============== CASE DATA ===============
## including new case data download here
## https://coronavirus.data.gov.uk/#
## this needs doing manually and coronavirus-cases_latest.csv placing in sitrep folder


fn <- glue(here::here('data'))+ '/UKD_' + td + '.Rdata'
if(file.exists(fn)){
  load(fn)
} else {
  UKD <- fread(here::here("../sitrep/coronavirus-cases_latest.csv"))
  save(UKD,file=fn)
}

load(here::here('data/pnts.Rdata'))     #this needs seeding
lcl <- 'Sheffield'                      #locale
lkd <- '23/03/2020'                     #lockdown date
lkdsig <- lkd #NEW '28/03/2020'                  #visible in data
(lkdnpt <- which(pnts$date==dmy(lkdsig)))
ndys <- length(UKD[,unique(`Specimen date`)])-lkdnpt #data fitted to most recent N days=days since lkdnpt
cat(ndys,file=here::here('data/ndys.txt'))


## cases and growth
## totals etc
UKD[,date:=`Specimen date`]
UKD[,type:=`Area type`]
UKD[,UTLA:=`Area name`]
UKD[,confirminc:=`Daily lab-confirmed cases`]
UKD[,confirm:=`Cumulative lab-confirmed cases`]

UKT <- UKD[type=='Country',.(cases=sum(confirm)),by=date]; UKT[,Population:='UK']
UKS <- UKD[UTLA==lcl,.(cases=sum(confirm)),by=date]; UKS[,Population:=lcl]
UKB <- rbind(UKT,UKS)


## doubling calculation NOTE decreasingly relevant
UKB[,date:=ymd(date)]
UKB[,dys:=date-min(date)]
growth <- UKB[date<=dmy(lkdsig),#date-days(ndys), NEW
{mod <- lm(log(1+cases)~dys);
  list(slope=coef(mod)['dys'],
       slope.lo=confint(mod,'dys',level=0.95)[1],
       slope.hi=confint(mod,'dys',level=0.95)[2])},
  by=Population]

growth[,doubling:=log(2)/slope];
growth[,doubling.lo:=log(2)/slope.hi];
growth[,doubling.hi:=log(2)/slope.lo]
growth <- merge(growth,UKB[,.(mxc=max(cases)),by=Population],by='Population')
growth[,loc:=mean(log(mxc))]
growth[,loc:=loc/2 + log(mxc)/2]
growth[,loc:=exp(loc)]                  #location for plot
growth[,loc:=30]                  #location for plot NEW
growth[,xloc:=as.Date(dmy(lkdsig) - days(14))]
growth[,txt:=format(round(doubling, 1), nsmall = 1)]
growth[,txt:=paste0('x2 every ',txt,' days\n(',
                    format(round(doubling.lo, 1), nsmall = 1),' to ',
                    format(round(doubling.hi, 1), nsmall = 1),')')]


## plot annotation locations & filename
xd <- dmy(lkd)+days(1)                  #marker for lockdown
xd <- format(xd,"%d/%m/%Y")
xd2 <- dmy(lkd)+days(2)                  #marker
xd2 <- format(xd2,"%d/%m/%Y")
pnm <- glue(here::here('plots')) + '/S1_' + td + '.pdf'       #

simp <- UKB[Population=='Sheffield']
simp <- simp[order(date)]
simp[,dcase:=c(0,diff(cases))]

GP1a <- ggplot(simp,
       aes(date,dcase,group=Population,col=Population)) +
  geom_line() +
  geom_point() +
  scale_y_log10(label=comma,breaks=log_breaks(n=10)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(0.25, 0.9),legend.direction='horizontal') + 
  annotate('rect',xmin=dmy(lkd),xmax=max(UKB$date),ymin=0,ymax=Inf,alpha=0.1)+
  annotation_logticks(sides='lr') +
  annotate('text',x=dmy(xd),y=1,label='lockdown') +
  geom_text(data=growth,aes(x=xloc,y=loc,label=txt,col=Population),
            show.legend = FALSE) +
  geom_vline(xintercept = dmy(lkdsig),col=2,lty=2)+ #ymd(UKB$date[lkdnpt+1])
  ggtitle('A) Growth rates in Sheffield and UK') +
  xlab('Date') + ylab('Cumulative confirmed COVID-19 cases (log scale)')
GP1a



## plot
GP1 <- ggplot(UKB,
       aes(date,cases,group=Population,col=Population)) +
  geom_line() +
  geom_point() +
  scale_y_log10(label=comma,breaks=log_breaks(n=10)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(0.25, 0.9),legend.direction='horizontal') + 
  annotate('rect',xmin=dmy(lkd),xmax=max(UKB$date),ymin=0,ymax=Inf,alpha=0.1)+
  annotation_logticks(sides='lr') +
  annotate('text',x=dmy(xd),y=1,label='lockdown') +
  geom_text(data=growth,aes(x=xloc,y=2*loc,label=txt,col=Population),
            show.legend = FALSE) +
  geom_vline(xintercept = dmy(lkdsig),col=2,lty=2)+ #ymd(UKB$date[lkdnpt+1])
  ggtitle('A) Growth rates in Sheffield and UK') +
  xlab('Date') + ylab('Cumulative confirmed COVID-19 cases (log scale)')
GP1

ggsave(GP1a,file=pnm,w=7,h=5)
ggsave(GP1a,file=here::here('figs/Trends1a.pdf'),w=7,h=5)
ggsave(GP1,file=here::here('figs/Trends.pdf'),w=7,h=5)


## make inference targets
UKS <- UKD[UTLA==lcl &  & `Area type`=='Upper tier local authority',
           .(cases=confirminc,date=ymd(date))]
UKS <- unique(UKS)
UKS[,dta:=TRUE]
UKS[,Population:='Sheffield']
dts <- min(UKS$date)
dts <- seq(from=max(as.Date(UKS$date)) + days(1),by='day',length.out = 100)
dts <- ymd(dts)
SF <- data.table(date=dts,#format(dts,"%d/%m/%Y"),
                 Population='Sheffield',
                 cases=NA,dta=FALSE)
SF <- rbind(UKS,SF,fill=TRUE)
SF[,dys:=(date)-min((date))]
ldy <- SF[!is.na(cases),max(dys)]
SF[,ddys:=dys-ldy]
ldt <- SF[ddys==0,cases]
who1 <- SF[,which(ddys==0)]              #last one
SF$cases <- as.numeric(SF$cases)

save(SF,file=here::here('data/SF.Rdata'))

## case targets
tmp <- UKD[UTLA=='Sheffield' & `Area type`=='Upper tier local authority']
tmp[,date:=ymd(date)]
pnts <- UKS[dta==TRUE]
pnts[,c('quantity','growth','grp'):=list('cases',NA,1)]
pnts[,date:=(date)]
pnts[,value:=cases]
pnts <- merge(pnts,tmp[,.(date,confirm)],by='date')

save(pnts,file=here::here('data/pnts.Rdata'))


 
## =============== SITREP DATA ================
## sitrep
D <- read_excel(here::here('../sitrep/DailySitrepMHman.xlsx'),skip = 2)
D <- as.data.table(D)
for (j in names(D))
  set(D,which(is.na(D[[j]])),j,0)

D[,`...7`:=NULL]
names(D)[1] <- 'date'
D[,date:=dmy(date)]

DM <- melt(D,id='date')

save(DM,file=here::here('data/DM.Rdata'))

## make calibration target file from these
vrs <- DM[,unique(variable)]
vrs
dss <- DM[variable==vrs[17],.(deaths=(value)),by=date] #NEW
tgts <- dss
nhsp <- DM[variable%in%vrs[12:13],.(newhosp=sum(value)),by=date] #NEW
tgts <- merge(tgts,nhsp,by='date',all.x=TRUE)
inhsp <- DM[variable%in%vrs[1:4],.(inhosp=sum(value)),by=date] #NOTE need to change
disch <- DM[variable==vrs[14],.(disch=(value)),by=date]        #NEW
tgts <- merge(tgts,disch,by='date',all.x=TRUE)
tgts <- merge(pnts[,.(date,cases)],tgts,by='date',all=TRUE)
tmp1 <- tgts[date>=max(date)-days(7)]    #keep NA in case of mismatch on last week
tgts <- tgts[date<max(date)-days(7)]
for (j in names(tgts))
  set(tgts,which(is.na(tgts[[j]])),j,0)
tgts <- rbind(tgts,tmp1)

save(tgts,file=here::here('data/tgts.Rdata'))

vrs

## lab testing
GPno <- ggplot(DM[variable %in% vrs[18:19]], #NEW
       aes(date,value,
           fill=factor(variable, levels=vrs[19:18]))) +
  geom_bar(stat='identity') +
  theme(legend.title = element_blank(),legend.position = 'top') + 
  ylab('Daily COVID lab tests') + xlab('Date') + xlim(c(dmy('20/03/2020'),NA))
GPno

GPno2 <- ggplot(DM[variable %in% vrs[18:19]], #NEW
       aes(date,value,
           fill=factor(variable, levels=vrs[19:18]))) +
  geom_bar(stat='identity',position='fill') +
  theme(legend.title = element_blank(),legend.position = 'top') + 
  ylab('Proportion of COVID lab tests') + xlab('Date') + xlim(c(dmy('20/03/2020'),NA))
GPno2


gb <- ggarrange(GPno,GPno2,common.legend = TRUE)
gb

ggexport(gb,filename=here::here('figs/LabTests.pdf'),height=5,width=7)

vrs                                     #NEW
vrs[14]                                 #disch
vrs[17]                                 #deaths
vrs[1:4]                                #beds
vrs[c(12,13)]                           #inflow

beds <- DM[variable %in% vrs[1:4],.(beds=sum(value)),by=date]
inout <- DM[variable %in% vrs[c(14,17)],.(outflow=sum(value)),by=date]
flw <- DM[variable %in% vrs[c(12,13)],.(inflow=sum(value)),by=date]
inout <- merge(inout,flw,by='date')
beds <- merge(beds,inout,by='date')
beds[,netflow:=inflow-outflow]
beds <- beds[order(date)]
beds[,netflow1:=c(0,rev(rev(netflow)[-1]))]
beds

## use net flow to update beds after:
bdy <- '2020-04-26'
initbed <- beds[date==bdy,beds]
beds[date>bdy,beds:=initbed+cumsum(netflow1)]
beds[date>bdy]

out <- dcast(DM[variable %in% vrs[c(14,17)]], #NEW
             date ~ variable,value.var = 'value')
out <- merge(out,beds,by='date')

out[,death.rate:=`Number of Patients Deceased with COVID 19 in Last 24HRS`/beds]
out[,discharge.rate:=`Number of COVID 19 Patients Discharges past 24HRS`/beds]
outm <- melt(out[,.(date,death.rate,discharge.rate)],id='date')

save(out,file=here::here('data/out.Rdata'))

ggplot(outm,aes(date,value,col=variable)) +
  geom_point() +  facet_wrap(~variable) +
  xlab('Date') + ylab('Raw daily rates') +
  geom_vline(xintercept = max(outm$date)-days(7),col=2,lty=2) +
  theme(legend.position = 'none')
## NOTE discharge rate clearly changin

ggsave('figs/OutRates.pdf',w=7,h=5)

## out rates
outr <- outm[,.(value=sum(value)),by=date]
stay <- outr[date>max(date)-days(7),.(mid=mean(1/value),sd=sd(1/value))]
stay <- as.list(stay)

## out rate
outr[is.finite(value),as.integer(date)]
md <- loess( value ~ as.integer(date), outr[is.finite(value)] )
outr[,v:=predict(md,newdata = as.integer(outr$date))]

GP <- ggplot(outr,aes(date,value)) + geom_point() + geom_line(aes(date,v)) + xlab('Date') + ylab('Rate of exiting hospital per day (death & discharge)')

ggsave(GP,file=here::here('figs/outrate.pdf'))



## stay fatality rates
n1 <- out[,sum(`Number of Patients Deceased with COVID 19 in Last 24HRS`)]
n2 <- out[,sum(`Number of COVID 19 Patients Discharges past 24HRS`)]
sfr <- list(sfr=n1/(n1+n2),n1=n1,n2=n2)

hosparms <- c(stay,sfr,outr=outr[!is.na(v)])
save(hosparms,file=here::here('data/parm.hosparms.Rdata'))
cat(round(stay$mid,d=1),file=here::here('data/meanstay.txt'))

## bed proportions
B <- DM[variable %in% vrs[1:4]]
B[,beds:=sum(value),by=date]

bpd1 <- ymd('2020-04-07')
bpd2 <- ymd('2020-04-21')

ggplot(B,aes(date,1e2*value/beds,col=variable)) + geom_line() + geom_point() +
  guides(col=guide_legend(nrow=2))+
  theme(legend.position = 'top',legend.title = element_blank()) +
  ylab('Percentage of beds by type') + xlab('Date') +
  xlim(c(dmy('20/03/2020'),NA)) +
  geom_vline(xintercept = bpd1,col=2,lty=2) +
  geom_vline(xintercept = bpd2,col=2,lty=2) 

ggsave(here::here('figs/BedProps.pdf'),w=7,h=5)

bedprops <- B[date>bpd1 & date<bpd2,.(prop=sum(value)/sum(beds)),by=variable]

save(bedprops,file=here::here('data/parm.bedprops.Rdata'))
fwrite(bedprops,file=here::here('data/parm.bedprops.csv'))

## O2 use
vrs
V <- DM[variable %in% vrs[9:11]]        #NEW
V <- merge(V,unique(B[,.(date,beds)]),by='date')

opd1 <- ymd('2020-04-09')
opd2 <- ymd('2020-04-25')


ggplot(V,aes(date,1e2*value/beds,col=variable)) + geom_line() + geom_point() +
  guides(col=guide_legend(nrow=3))+
  theme(legend.position = 'top',legend.title = element_blank()) +
  ylab('Percentage of breathing assistance by type') + xlab('Date')+
  xlim(c(dmy('20/03/2020'),NA)) +
  geom_vline(xintercept = opd1,col=2,lty=2) +
  geom_vline(xintercept = opd2,col=2,lty=2) 

ggsave(here::here('figs/Oprops.pdf'),w=7,h=5)


o2props <- V[date>opd1 & date <opd2,.(prop=sum(value)/sum(beds)),by=variable]

save(o2props,file=here::here('data/parm.o2props.Rdata'))
fwrite(o2props,file=here::here('data/parm.o2props.csv'))


## vrs
## delays cases, hosp, deaths
admns <- DM[variable %in% vrs[c(12,13)],.(mid=sum(value)),by=date] #NEW
admns[,quantity:='hosp']
deaths <- DM[variable==vrs[17],.(mid=value),by=date] #NEW
deaths[,quantity:='deaths']

## load(here::here('data/pnts.Rdata'))

inc <- rbind(admns,deaths)
inc <- rbind(pnts[,.(date,mid=cases,quantity)],inc)
inc[,cinc:=cumsum(mid),by=quantity]

## winc <- dcast(inc[cinc<100],date~quantity,value.var = 'cinc')
## winc[,ndys:=yday(date)]
## winc[,ndys:=ndys-ndys[1]]

gc <- ggplot(inc,aes(date,cinc,col=quantity)) +
  geom_line() + geom_point() +
  scale_y_log10() + ylab('Cumulative numbers') + xlab('Date')

cfs <- coef(lm(data=inc[cinc<100],cinc~date + quantity))

(dC2D <- cfs[3]/-cfs[2])
(dC2H <- cfs[4]/-cfs[2])

inc[,date2:=date]
inc[quantity=='deaths',date2:=date - days(round(dC2D))]
inc[quantity=='hosp',date2:=date - days(round(dC2H))]

gc2 <- ggplot(inc,aes(date2,cinc,col=quantity)) +
  geom_line() + geom_point() +
  ## geom_smooth(method='lm')+
  scale_y_log10() + ylab('Cumulative numbers') + xlab('Shifted date')

gcb <- ggarrange(gc,gc2,common.legend = TRUE)
gcb

ggexport(gcb,filename=here::here('figs/Delays.pdf'),width=7,height=5)

delays <- list(dC2D=dC2D,dC2H=dC2H)
save(delays,file=here::here('data/parm.delays.Rdata'))
cat(round(delays$dC2D,d=1),file=here::here('data/dC2D.txt'))
cat(round(delays$dC2H,d=1),file=here::here('data/dC2H.txt'))
