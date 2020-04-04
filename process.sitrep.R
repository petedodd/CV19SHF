## process sitrep
library(here)
library(readxl)
library(data.table)
library(lubridate)
library(ggpubr)

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

vrs <- DM[,unique(variable)]

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

## stay fatality rates
n1 <- out[,sum(`Number of Patients Deceased with COVID 19 in Last 24HRS`)]
n2 <- out[,sum(`Number of COVID 19 Patients Discharges past 24HRS`)]
sfr <- list(sfr=n1/(n1+n2),n1=n1,n2=n2)

hosparms <- c(stay,sfr)

save(hosparms,file=here::here('data/parm.hosparms.Rdata'))
cat(stay$mid,file=here::here('data/meanstay.txt'))

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

load(here::here('data/pnts.Rdata'))

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
cat(delays$dC2D,file=here::here('data/dC2D.txt'))
cat(delays$dC2H,file=here::here('data/dC2H.txt'))
