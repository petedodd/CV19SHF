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

lcl <- 'Sheffield'                      #locale
lkd <- '23/03/2020'                     #lockdown date
lkdsig <- '28/03/2020'                  #visible in data
lkdnpt <- 21                #before intervention signal (or consider 16=actual lockdown)
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
growth <- UKB[date>=dmy(lkdsig),#date-days(ndys),
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
  geom_text(data=growth,aes(x=dmy(xd2),y=loc,label=txt,col=Population),
            show.legend = FALSE) +
  geom_vline(xintercept = dmy(lkdsig),col=2,lty=2)+ #ymd(UKB$date[lkdnpt+1])
  ggtitle('A) Growth rates in Sheffield and UK') +
  xlab('Date') + ylab('Cumulative confirmed COVID-19 cases (log scale)')
GP1

ggsave(GP1,file=pnm,w=7,h=5)
ggsave(GP1,file=here::here('figs/Trends.pdf'),w=7,h=5)


## make inference targets
UKS <- UKD[UTLA==lcl,.(cases=confirminc,date=ymd(date))]
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
tmp <- UKD[UTLA=='Sheffield']
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
dss <- DM[variable==vrs[13],.(deaths=(value)),by=date]
tgts <- dss
nhsp <- DM[variable%in%vrs[8:9],.(newhosp=sum(value)),by=date]
tgts <- merge(tgts,nhsp,by='date',all.x=TRUE)
inhsp <- DM[variable%in%vrs[1:4],.(inhosp=sum(value)),by=date]
disch <- DM[variable==vrs[10],.(disch=(value)),by=date]
tgts <- merge(tgts,disch,by='date',all.x=TRUE)
tgts <- merge(pnts[,.(date,cases)],tgts,by='date',all=TRUE)
tmp1 <- tgts[date>=max(date)-days(7)]    #keep NA in case of mismatch on last week
tgts <- tgts[date<max(date)-days(7)]
for (j in names(tgts))
  set(tgts,which(is.na(tgts[[j]])),j,0)
tgts <- rbind(tgts,tmp1)

save(tgts,file=here::here('data/tgts.Rdata'))

## lab testing
GPno <- ggplot(DM[variable %in% vrs[14:15]],
       aes(date,value,
           fill=factor(variable, levels=vrs[15:14]))) +
  geom_bar(stat='identity') +
  theme(legend.title = element_blank(),legend.position = 'top') + 
  ylab('Daily COVID lab tests') + xlab('Date') + xlim(c(dmy('20/03/2020'),NA))
GPno

GPno2 <- ggplot(DM[variable %in% vrs[14:15]],
       aes(date,value,
           fill=factor(variable, levels=vrs[15:14]))) +
  geom_bar(stat='identity',position='fill') +
  theme(legend.title = element_blank(),legend.position = 'top') + 
  ylab('Proportion of COVID lab tests') + xlab('Date') + xlim(c(dmy('20/03/2020'),NA))
GPno2


gb <- ggarrange(GPno,GPno2,common.legend = TRUE)
gb

ggexport(gb,filename=here::here('figs/LabTests.pdf'),height=5,width=7)


vrs[10]                                 #disch
vrs[13]                                 #deaths
vrs[1:4]                                #beds

beds <- DM[variable %in% vrs[1:4],.(beds=sum(value)),by=date]
out <- dcast(DM[variable %in% vrs[c(10,13)]],date ~ variable,value.var = 'value')
out <- merge(out,beds,by='date')

out[,death.rate:=`Number of Patients Deceased with COVID 19 in Last 24HRS`/beds]
out[,discharge.rate:=`Number of COVID 19 Patients Discharges past 24HRS`/beds]
outm <- melt(out[,.(date,death.rate,discharge.rate)],id='date')

ggplot(outm,aes(date,value,col=variable)) +
  geom_point() +  facet_wrap(~variable) +
  xlab('Date') + ylab('Raw daily rates') +
  geom_vline(xintercept = max(outm$date)-days(7),col=2,lty=2) +
  theme(legend.position = 'none')

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

ggplot(B,aes(date,1e2*value/beds,col=variable)) + geom_line() + geom_point() +
  guides(col=guide_legend(nrow=2))+
  theme(legend.position = 'top',legend.title = element_blank()) +
  ylab('Percentage of beds by type') + xlab('Date') +
  xlim(c(dmy('20/03/2020'),NA)) +
  geom_vline(xintercept = max(outm$date)-days(7),col=2,lty=2) 

ggsave(here::here('figs/BedProps.pdf'),w=7,h=5)

bedprops <- B[date>max(date)-days(7),.(prop=sum(value)/sum(beds)),by=variable]

save(bedprops,file=here::here('data/parm.bedprops.Rdata'))
fwrite(bedprops,file=here::here('data/parm.bedprops.csv'))

## O2 use
V <- DM[variable %in% vrs[5:7]]
V <- merge(V,unique(B[,.(date,beds)]),by='date')


ggplot(V,aes(date,1e2*value/beds,col=variable)) + geom_line() + geom_point() +
  guides(col=guide_legend(nrow=3))+
  theme(legend.position = 'top',legend.title = element_blank()) +
  ylab('Percentage of breathing assistance by type') + xlab('Date')+
  xlim(c(dmy('20/03/2020'),NA)) +
  geom_vline(xintercept = max(outm$date)-days(7),col=2,lty=2) 

ggsave(here::here('figs/Oprops.pdf'),w=7,h=5)


o2props <- V[date>max(date)-days(7),.(prop=sum(value)/sum(beds)),by=variable]

save(o2props,file=here::here('data/parm.o2props.Rdata'))
fwrite(o2props,file=here::here('data/parm.o2props.csv'))


## delays cases, hosp, deaths
admns <- DM[variable %in% vrs[8:9],.(mid=sum(value)),by=date]
admns[,quantity:='hosp']
deaths <- DM[variable==vrs[13],.(mid=value),by=date]
deaths[,quantity:='deaths']

## load(here::here('data/pnts.Rdata'))

inc <- rbind(admns,deaths)
inc <- rbind(pnts[,.(date,mid=cases,quantity)],inc)
inc[,cinc:=cumsum(mid),by=quantity]

gc <- ggplot(inc,aes(date,cinc,col=quantity)) +
  geom_line() + geom_point() +
  ## geom_smooth(method='lm')+
  scale_y_log10() + ylab('Cumulative numbers') + xlab('Date')

cfs <- coef(lm(data=inc,cinc~date + quantity))

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
