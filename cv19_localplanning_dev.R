## A simple analysis for local COVD-19 capacity forecasting
## NB does not predict intervention effect, except via indications in data
## 28 March 2020
## p.j.dodd@sheffield.ac.uk

## --- load Libraries
library(tidyr)
library(ggplot2)
library(ggpubr)
library(data.table)
library(here)
library(glue)
library(lubridate)
library(scales)
library(odin)
library(MASS)
library(curl)

## --- set-up & utilities
## make sub directories
if(!file.exists(here::here('data'))) dir.create(here::here('data'),showWarnings=FALSE)
if(!file.exists(here::here('plots'))) dir.create(here::here('plots'),showWarnings=FALSE)
getdate <- function()glue(gsub('-','_',Sys.Date()))
td <- getdate()                         #date stamp

## for color consistency
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

shfc <- gg_color_hue(2)[1]              #sheffield color for consistency

## get data
## https://github.com/emmadoughty/Daily_COVID-19

fn <- glue(here::here('data'))+ '/UKC_' + td + '.Rdata'
if(file.exists(fn)){
  load(fn)
} else {
  UKC <- fread("https://raw.githubusercontent.com/emmadoughty/Daily_COVID-19/master/Data/cases_by_utla.csv")
  save(UKC,file=fn)
}


## UKD <- fread("https://raw.githubusercontent.com/emmadoughty/Daily_COVID-19/master/Data/deaths_by_area.csv")

lcl <- 'Sheffield'                      #locale
lkd <- '23/03/2020'                     #lockdown date
lkdnpt <- 21                #before intervention signal (or consider 16=actual lockdown)
ndys <- length(UKC[,unique(date)])-lkdnpt #data fitted to most recent N days=days since lkdnpt


## parameters etc for later
## demography etc from Chris Gibbon's work book
ages <- c('0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+')
shfdemo <- c(0.107886922,0.115624544,0.188336768,0.139600136,0.11854238,0.122764966,
             0.090433456,0.072074455,0.044736372)
sympto <- c(66.00,66.00,66.00,66.00,66.00,66.00,66.00,66.00,66.00)
symptohosp <- c(0.10,0.30,1.20,3.20,4.90,10.20,16.60,24.30,27.30)
hospcc <- c(5.00,5.00,5.00,5.00,6.30,12.20,27.40,43.20,70.90)
IFR <- c(0.00,0.01,0.00,0.08,0.15,0.60,2.20,5.10,9.30)

parms <- data.table(ages=ages,demo=shfdemo,sympto=sympto,symptohosp=symptohosp,
                    hospcc=hospcc,IFR=IFR)
parms[,proph:=sympto*symptohosp/1e4]    #proportion hospitalized
parms[,propcc:=sympto*symptohosp*hospcc/1e6] #proportion crit care
parms[,propd:=IFR/1e2]                       #proportion who will die

## delay and parameter delays
D2D <- days(20)                        #delay to death
D2CC <- D2H <- days(10)                 #delay to hosp
mn.hosp.stay <- 7                               #mean hosp stay
mn.cc.stay <- 10                                #mean CC stay
gt <- 6.5                               #generation time in days

## mean parms for Sheffield pop
proph <- parms[,weighted.mean(proph,w=demo)]
propcc <- parms[,weighted.mean(propcc,w=demo)]
propd <- parms[,weighted.mean(propd,w=demo)]


## mm <- lm(data=UKB[dmy(date)>=max(dmy(date))-days(ndys)],log(1+cases) ~dys)

## totals etc
UKT <- UKC[type=='Country',.(cases=sum(confirm)),by=date]; UKT[,Population:='UK']
UKS <- UKC[UTLA==lcl,.(cases=sum(confirm)),by=date]; UKS[,Population:=lcl]
UKB <- rbind(UKT,UKS)


## doubling calculation
UKB[,date:=dmy(date)]
UKB[,dys:=date-min(date)]
growth <- UKB[date>=date-days(ndys),
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
pnm <- glue(here::here('plots')) + '/P1_' + td + '.pdf'       #


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
  ## annotate('text',x=dmy(xd),y=0.25*growth$loc[1],col=2,hjust=1,
  ##          label='Lockdown effects may be yet to appear') +
  geom_text(data=growth,aes(x=dmy(xd2),y=loc,label=txt,col=Population),
            show.legend = FALSE) +
  ggtitle('A) Growth rates in Sheffield and UK') +
  xlab('Date') + ylab('Daily confirmed COVID-19 cases (log scale)')
GP1

ggsave(GP1,file=pnm,w=7,h=5)

## difference with start added & -ve -> 0
df1 <- function(x){
  x <- c(x[1],diff(x))
  x[x<0] <- 0
  x
}

## make data for sheffield prediction
UKB[,cases:=df1(cases),by=Population]   #introduce diff
UKB <- UKB[dmy(UKB$date)!=min(dmy(UKB$date))]

UKS <- UKB[Population=='Sheffield']
UKS[,dta:=TRUE]
dts <- min(UKS$date)
dts <- seq(from=max(dmy(UKS$date)) + days(1),by='day',length.out = 100)
SF <- data.table(date=format(dts,"%d/%m/%Y"),Population='Sheffield',
                 cases=NA,dta=FALSE)
SF <- rbind(UKS,SF,fill=TRUE)
SF[,dys:=dmy(date)-min(dmy(date))]
## sheffield pop
shfpop <- 616210
## md <- lm(data=SF[as.Date(date)>=lkd & !is.na(cases)],log(1+cases) ~ dys)
md <- lm(data=SF[dmy(date)>=ymd(td)-days(ndys) & !is.na(cases)],
         log(1+cases) ~ dys)
ldy <- SF[!is.na(cases),max(dys)]
SF[,ddys:=dys-ldy]
ldt <- SF[ddys==0,cases]
who1 <- SF[,which(ddys==0)]              #last one
cfs <- c(coef(md)['dys'],confint(md,'dys',level=0.95)[1:2])
SF$cases <- as.numeric(SF$cases)


## point data
pnts <- UKS[dta==TRUE]
pnts[,c('quantity','growth','grp'):=list('cases',NA,1)]
pnts[,date:=dmy(date)]
pnts[,value:=cases]

tz <- seq(from=0,to=100,by=1)
css <- pnts[,(cases)]                   #t = dys in this data
Rinit <- pnts[,sum(cases)]

UR <- 1                               #underreporting

## interventions assessment
lat <- 5.1                              #latent
pinf <- 4.6                              #infectious

## SEIR model
seiri <- odin::odin({
  initial(S) <- N-I0 - Rinit
  initial(E) <- I0/2
  initial(I) <- I0/2
  initial(R) <- Rinit
  deriv(S) <- -beta*S*I/N
  deriv(E) <- beta*S*I/N - E*nu
  deriv(I) <- E*nu - I*nu2
  deriv(R) <- +I*nu2
  output(incidence) <- E*nu
  nu <- 1/lat
  nu2 <- 1/pinf
  beta <- R0*nu2*z
  lat <- user(5.1)                       #latent time
  pinf <- user(4.6)                      #infectious period
  I0 <- user(361)                        #initial state
  Rinit <- user(5)                      #initial recovered
  R0 <- user(2.6)
  N <- user(616210)
  z <- interpolate(tt, y)
  tt[] <- user()
  y[] <- user()
  dim(tt) <- user()
  dim(y) <- length(tt)
},target='r')

## likelihood

HR <- 0.5                               #hazard ratio of effect
y <- tz
y[1:lkdnpt] <- 1
y[(1+lkdnpt):length(y)] <- HR

## ## test
## modi <- siri(I0=exp(res2$par[1]),R0=exp(res2$par[2]),Rinit = Rinit/UR,y=y,tt=tz)
## out3 <- modi$run(tz)
## plot(out3)
## par(mfrow=c(1,1))


seiriLL <- function(x){
  y[1:lkdnpt] <- 1
  y[(1+lkdnpt):length(y)] <- exp(x[3])
  mod <- seiri(I0=exp(x[1]),R0=exp(x[2]),Rinit = Rinit/UR,y=y,tt=tz)
  out <- mod$run(tz)
  ecss <- out[1:length(css),'incidence']
  sum(dpois(css,UR*ecss,log=TRUE)) #last ndys data
}

## seiriLL(c(0,1,-0.1))                    #test
## pnts

## compare sheffield
UR <- 0.075 #UK reporting: https://cmmid.github.io/topics/covid19/severity/global_cfr_estimates.html
resi <- optim(par=c(0,1,0),fn=seiriLL,control = list(fnscale=-1),hessian=TRUE) #ML

y[1:lkdnpt] <- 1
y[(1+lkdnpt):length(y)] <- exp(resi$par[3])
mod <- seiri(I0=exp(resi$par[1]),R0=exp(resi$par[2]),Rinit = Rinit/UR,y=y,tt=tz)
outi <- mod$run(tz)


## plot(outi)
## par(mfrow=c(1,1))

## plot(1:length(css),css,col=c(rep(2,lkdnpt),rep(1,(length(css)-lkdnpt))))

## check fit
ylmz <- c(1,max(css,UR*outi[1:length(css),'incidence']))

## plot(1:length(css),UR*outi[1:length(css),'incidence'],type='l',log='y',ylim = ylmz,
##      xlab='days',ylab='log daily confirmed cases')
## points(1:length(css),css,col=c(rep(2,lkdnpt),rep(1,(length(css)-lkdnpt))))


## real scale
pdf(here::here('plots/FitExplain.pdf'))

plot(1:length(css),UR*outi[1:length(css),'incidence'],type='l',
     ylim=c(0,max(css,UR*outi[1:length(css),'incidence'])),
     xlab='days',ylab='log daily confirmed cases')
points(1:length(css),css,col=c(rep(2,lkdnpt),rep(1,(length(css)-lkdnpt))))

text(10,40,'Red points fix R0 & initial state',col=2)
text(10,80,'Black points fix intervention effect')

dev.off()

## extrapolation
plot(tz,UR*outi[,'incidence'],type='l')
points(1:length(css),css,col=c(rep(2,lkdnpt),rep(1,(length(css)-lkdnpt))))


pntsd <- data.table(id=1,t=1:nrow(pnts),notes=pnts$cases)
LA <- data.table(notes = UR*outi[,'incidence'],trueincidence=outi[,'incidence'],t=tz)

PP <- ggplot(LA,aes(x=t,y=notes)) +
  geom_line()  +
  geom_point(data=pntsd) +
  geom_line(aes(t,trueincidence),col=2)
PP
LA[,max(trueincidence)]
LA[,sum(trueincidence)]

PP + scale_y_log10()


## ## checks and uncertainty
## ## this is rather an under-estimate
## ## uncertainty over parameter estimates
## Sig <- solve(-resi$hessian)
## PMZ <- mvrnorm(200,resi$par,Sigma=Sig)
## LA <- list()
## for(i in 1:nrow(PMZ)){
##   y[1:lkdnpt] <- 1; y[(1+lkdnpt):length(y)] <- exp(resi$par[3])
##   mod <- seiri(I0=exp(resi$par[1]),R0=exp(resi$par[2]),Rinit = Rinit/UR,y=y,tt=tz)
##   outi <- mod$run(tz)
##   LA[[i]] <- as.data.table(mod$run(tz))
## }

## ggplot(LA,aes(x=t,y=incidence*UR,group=id)) +
##   geom_line()  +
##   geom_point(data=pntsd,aes(t,notes)) +
##   scale_y_log10()


## LAM <- LA[,.(mid=median(incidence),
##             lo=quantile(incidence,.025),
##             hi=quantile(incidence,.975)),
##          by=t]

## LA <- rbindlist(LA)
## LA[,id:=rep(1:length(tz),each=nrow(PMZ))]
## LA

LA
SF
ou3 <- merge(SF,LA,by.x='dys',by.y='t')

UKS <- ou3[,.(date,
              cases=notes,
              truecases=trueincidence,
              dys)]


## calculate & delays
UKS[,hosp:=proph*truecases];
UKS[,ccadm:=propcc*truecases]; 
UKS[,deaths:=propd*truecases]; 

## reformat
SM <- melt(UKS[,.(date=dmy(date),
                  cases,truecases,hosp,ccadm,deaths
                  )],id.vars = 'date')

names(SM)[2] <- 'quantity'

## introduce delays
SM[quantity=='deaths',date:=date + D2D]
SM[quantity=='hosp',date:=date + D2H]
SM[quantity=='ccadm',date:=date + D2CC]

## case projections on real scale
lkdl <- dmy(lkd)
xdl <- dmy(xd)
pnm <- glue(here::here('plots')) + '/S2_' + td + '.pdf'       #

## incidence plot
GP3 <- ggplot(SM[date<=dmy('15/07/2020')],
              aes(date,value,col=quantity)) +
  geom_point(data=pnts) +
  geom_line()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  xlab('Date') + ylab('Daily incidence of quantity (log scale)') +
  scale_y_log10(label=comma,breaks=log_breaks(n=10),limits=c(1,NA))  +
  scale_color_manual(breaks=c('truecases','cases','hosp','ccadm','deaths'),
                     labels=c('true incidence','confirmed cases','hospital admissions',
                              'critical care admissions','deaths'),
                     values=cbbPalette[1:5])+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(0.85, 0.5),legend.direction='vertical') + 
  annotation_logticks(sides='lr') +
  geom_vline(xintercept = ymd(pnts$date[lkdnpt]),col=2,lty=2)+
  annotate('text',x=ymd(pnts$date[lkdnpt])-days(5),y=100,col=2,
           label=paste0('R0 = ',round(exp(resi$par[2]),digits=1)))+
  annotate('text',x=ymd(pnts$date[lkdnpt])+days(10),y=140,col=2,
           label=paste0('transmission reduction = ',
                        round((1-exp(resi$par[3]))*1e2,digits=0),"%"))+
  ggtitle('B) Projected daily numbers for Sheffield')

GP3

## parameter table
GT <- ggtexttable(parms[,.(ages,
                           `% of people`=format(round(1e2*demo,1),nsmall=1),
                           `% symptomatic`=sympto,
                           `% of symptomatic\nhospitalised`=symptohosp,
                           `% of hospitalised\ncritical care`=hospcc,
                           `IFR (%)`=IFR)],
                  rows = NULL,
                  theme = ttheme("classic",base_size=5))

## assumptions text
TXT <- glue('ASSUMPTIONS:')
TXT <- TXT + '\nAssumptions on outcomes left'
TXT <- TXT + '\nDemography used to compute average outcomes\n'
TXT <- TXT + '\nMean delays case:(hospitalisation, critical care, death ) = \n('
TXT <- TXT + as.duration(D2H)/ddays(1) + ', ' + as.duration(D2CC)/ddays(1)
TXT <- TXT + ', ' + as.duration(D2D)/ddays(1) + ') days\n'
TXT <- TXT + '\nMean stays:(hospitalisation, critical care) = \n('
TXT <- TXT + mn.hosp.stay + ', ' + mn.cc.stay + ') days'
TXT <- TXT + '\nProjection of incidence via SEIR model:'
TXT <- TXT + '\n - change in beta on dashed line = ' + round(1e2*(1-exp(resi$par[3])))+ '% redn'
TXT <- TXT + '\n - 7.5% infections confirmed'
TXT <- TXT + '\n - (latent,infectious) periods = ( ' + lat + ',' + pinf + ') days'
GE <- ggplot(data.frame(x=c(0,1),y=c(0,1)),aes(x,y)) +
  geom_point(col=NA)+
  annotate('text',x=0,y=0.5,col=2,hjust = 0,label=TXT,size=3) +
  theme_void()
GE

GA <- ggarrange(GP3,ggarrange(GT,GE,nrow=1,ncol=2),nrow=2,ncol=1,heights=c(2,1))

## panel B on incidences
pnm <- glue(here::here('plots')) + '/S3_' + td + '.pdf'       #
ggexport(GA,filename = pnm)


## prevalence
ldd <- max(pnts$date)
nac <- function(x){x[is.na(x)] <- 0;cumsum(x)}
PM <- copy(SM)
PM <- PM[order(quantity,date)]
## calculate cumulatives
PM[quantity%in%c('truecases','deaths'),
   value:=nac(value),by=.(quantity)]

## prevalence as simple incidence x duration
PM[quantity=='hosp',value:=mn.hosp.stay*(value),by=.(quantity)]
PM[quantity=='ccadm',value:=mn.cc.stay*(value),by=.(quantity)]

## max labels for graph
MX <- PM[,.(value=max(value,na.rm=TRUE)),by=quantity]
hmax <- MX[quantity=='hosp',value]
dmax <- PM[ quantity=='hosp' & value==hmax,date]
dpk <- SM[quantity=='deaths',max(value)]
dpk <- format(signif(dpk,3),big.mark=',')
MX[,value:=format(signif(value,3),big.mark=',')]
MX <- MX[order(quantity)]
MX <- MX[quantity!='cases']
MX[,lbl:=c('cumulative cases','hospital beds\n(including critical care)',
           'critical care beds','cumulative deaths'
           )]
MX[,lbl:=paste0(value,' ',lbl)]
MX[quantity=='deaths',lbl:=paste0(lbl," (",dpk,"/day peak)")]
MX[,date:=rep(dmax,4)]
MX[,value:=rep(hmax,4) + 0e3]
MX[,value:=value * c(3^3,3^2,3^{1},3^0)]
## MX[,grp:='truecases1'];MX[,cases.over.confirmed:='1']
MX

dmid <- min(PM$date) + days(max(PM$date)-min(PM$date))/2
MX$date <- dmid

## prevalence plot
GP4 <- ggplot(PM[date<=dmy('15/07/2020') & value>0 & quantity!='cases'],
              aes(date,value,col=quantity)) +
  geom_line() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  xlab('Date') + ylab('Prevalence of quantity (log scale)') +
  scale_y_log10(label=comma,breaks=log_breaks(n=10),limits=c(1,NA)) +
  geom_text(data=MX,aes(label=lbl),show.legend = FALSE,hjust=0) + 
  scale_color_manual(breaks=c('truecases','hosp','ccadm','deaths'),
                     labels=c('cumulative cases',
                              'hospital beds',
                              'critical care beds','cumulative deaths'),
                     values=cbbPalette[1:4])+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(0.15, 0.7),legend.direction='vertical') +
  annotate('text',x=dmid,y=MX[,exp(mean(log(value)))],label='Max:',hjust=1.5,col=2)+
  annotation_logticks(sides='lr') +
  annotate('segment',x=dmax,xend=dmax,yend=1,y=hmax,lty=2,col=2)+
  annotate('text',x=dmax,y=1,label=paste0(dmax),hjust=1,col=2)+
  ggtitle('C) Projected prevalent or cumulative numbers for Sheffield') +
  expand_limits(x=PM[,min(date)])
GP4

pnm <- glue(here::here('plots')) + '/S4_' + td + '.pdf'       #
ggsave(GP4,filename = pnm,w=8,h=6)



All <- ggarrange(GP1,GA,GP4,ncol=1)
pnm <- glue(here::here('plots')) + '/All_SEIR' + td + '.pdf'       #
ggexport(All,filename=pnm)

All2 <- ggarrange(GP1,GA,GP4,ncol=1,nrow=3,heights = c(0.9,1.25,1))
pnm <- glue(here::here('plots')) + '/SAll_SEIR_' + td + '.pdf'       #
ggexport(All2,filename=pnm,width=7,height = 18)

## save out data too
pnm <- glue(here::here('data')) + '/SPM_' + td + '.Rdata'       #
save(PM,file = pnm)
pnm <- glue(here::here('data')) + '/SSM_' + td + '.Rdata'       #
save(SM,file = pnm)
