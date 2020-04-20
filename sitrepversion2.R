## new version to use additional private data
## further altered to use new PHE case data
## before running this, also run the sitrep updates
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
library(ggrepel)

## --- set-up & utilities
## make sub directories
if(!file.exists(here::here('data'))) dir.create(here::here('data'),showWarnings=FALSE)
if(!file.exists(here::here('plots'))) dir.create(here::here('plots'),showWarnings=FALSE)
if(!file.exists(here::here('figs'))) dir.create(here::here('figs'),showWarnings=FALSE)
getdate <- function()glue(gsub('-','_',Sys.Date()))
td <- getdate()                         #date stamp

## for color consistency
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

shfc <- gg_color_hue(2)[1]              #sheffield color for consistency

## sitrep  & case targets
load(here::here('data/DM.Rdata'))
vrs <- DM[,unique(variable)]
load(here::here('data/tgts.Rdata'))
load(here::here('data/pnts.Rdata'))

tgts
pnts

## sheffield pop
shfpop <- 616210

lcl <- 'Sheffield'                      #locale
lkd <- '23/03/2020'                     #lockdown date
lkdsig <- '28/03/2020'                  #visible in data
#before intervention signal (or consider 16=actual lockdown)
(lkdnpt <- which(pnts$date==dmy(lkdsig)))

## parameters etc for later
## sitrep parms
load(here::here('data/parm.delays.Rdata')) # from sitrep
load(here::here('data/parm.bedprops.Rdata')) #from sitrep
load(here::here('data/parm.hosparms.Rdata')) #from sitrep
load(here::here('data/parm.o2props.Rdata')) #from sitrep

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

## mean parms for Sheffield pop
proph <- parms[,weighted.mean(proph,w=demo)]
propcc <- parms[,weighted.mean(propcc,w=demo)]
propd <- parms[,weighted.mean(propd,w=demo)]

## ============== new ODE version
## partially to cope with censoring and non-eqm
seirid <- odin::odin({
  initial(S) <- N-I0 - Rinit
  initial(E) <- I0/2
  initial(I) <- I0/2
  initial(R) <- Rinit
  initial(dying) <- D0
  initial(sick) <- K0
  initial(hosp) <- 0
  deriv(S) <- -beta*S*I/N
  deriv(E) <- beta*S*I/N - E*nu
  deriv(I) <- E*nu - I*nu2
  deriv(R) <- +I*nu2
  deriv(sick) <- E*nu*HFR - sick/d2h    #waiting state for hospitalised fraction
  deriv(dying) <- E*nu*IFR - dying/d2d #copy state to capture delay to death
  deriv(hosp) <- sick/d2h - hosp/h2o      #prevalence in hosp
  output(incidence) <- E*nu             #linked to case detection
  output(deaths) <- dying/d2d
  output(admns) <- sick/d2h
  output(disch) <- hosp/h2o
  nu <- 1/lat
  nu2 <- 1/pinf
  beta <- R0*nu2*z
  lat <- user(5.1)                       #latent time
  pinf <- user(2.9)                      #infectious period
  I0 <- user(361)                        #initial state
  Rinit <- user(5)                      #initial recovered
  R0 <- user(2.6)
  N <- user(616210)
  d2d <- user(14)                         #delay case to death
  d2h <- user(10)                         #delay to hospn
  h2o <- user(14)                         #
  IFR <- user(1e-2)                       #IFR
  HFR <- user(1e-2)                       #hosp fraction
  D0 <- user(1)
  K0 <- user(1)
  z <- interpolate(tt, y)
  tt[] <- user()
  y[] <- user()
  dim(tt) <- user()
  dim(y) <- length(tt)
},target='r')

## ## version with dynamic stay time
## seirid2 <- odin::odin({
##   initial(S) <- N-I0 - Rinit
##   initial(E) <- I0/2
##   initial(I) <- I0/2
##   initial(R) <- Rinit
##   initial(dying) <- D0
##   initial(sick) <- K0
##   initial(hosp) <- 0
##   deriv(S) <- -beta*S*I/N
##   deriv(E) <- beta*S*I/N - E*nu
##   deriv(I) <- E*nu - I*nu2
##   deriv(R) <- +I*nu2
##   deriv(sick) <- E*nu*HFR - sick/d2h    #waiting state for hospitalised fraction
##   deriv(dying) <- E*nu*IFR - dying/d2d #copy state to capture delay to death
##   deriv(hosp) <- sick/d2h - hosp/h2o      #prevalence in hosp
##   output(incidence) <- E*nu             #linked to case detection
##   output(deaths) <- dying/d2d
##   output(admns) <- sick/d2h
##   output(disch) <- hosp/h2o
##   nu <- 1/lat
##   nu2 <- 1/pinf
##   beta <- R0*nu2*z
##   lat <- user(5.1)                       #latent time
##   pinf <- user(2.9)                      #infectious period
##   I0 <- user(361)                        #initial state
##   Rinit <- user(5)                      #initial recovered
##   R0 <- user(2.6)
##   N <- user(616210)
##   d2d <- user(14)                         #delay case to death
##   d2h <- user(10)                         #delay to hospn
##   IFR <- user(1e-2)                       #IFR
##   HFR <- user(1e-2)                       #hosp fraction
##   D0 <- user(1)
##   K0 <- user(1)
##   z <- interpolate(tt, y)
##   h2o <- interpolate(tt, y2o)
##   tt[] <- user()
##   y[] <- user()
##   y2o[] <- user()
##   dim(tt) <- user()
##   dim(y) <- length(tt)
##   dim(y2o) <- length(tt)
## },target='r')


## testing
tz <- seq(from=0,to=120,by=1)
HR <- 0.5                               #hazard ratio of effect
y <- tz
y[1:lkdnpt] <- 1
y[(1+lkdnpt):length(y)] <- HR
css <- pnts[,(cases)]                   #t = dys in this data
Rinit <- pnts[,sum(cases)]
y2o <- tz                               #stay dynamic
start <- which(tgts$date==hosparms$outr.date[1])
end <- which(tgts$date==rev(hosparms$outr.date)[1])
y2o[1:start] <- 1/(hosparms$outr.v[1]+5e-4)
y2o[end:length(tz)] <- 1/rev(hosparms$outr.v)[1]
y2o[start:end] <- 1/(hosparms$outr.v+5e-4)


## test
UR <- 1                                 #underreporting
mdi <- seirid(I0=exp(1),R0=4,Rinit = Rinit/UR,y=y,tt=tz,
               d2h=10,d2d=14,h2o=10,IFR=propd,HFR=proph)

## owi <- mdi$run(tz)
## plot(owi)
## par(mfrow=c(1,1))

css
tgts$cases

## run doer
dorun <- function(x){
  y[1:lkdnpt] <- 1
  y[(1+lkdnpt):length(y)] <- exp(x[3])
  md <- seirid(I0=exp(x[1]),R0=exp(x[2]),Rinit = Rinit/UR,y=y,tt=tz,
               d2h=delays$dC2H,d2d=delays$dC2D,h2o=hosparms$mid,
               IFR=propd,HFR=exp(x[4]))
  md$run(tz)
}

## mdi <- seirid2(I0=exp(1),R0=4,Rinit = Rinit/UR,y=y,tt=tz,y2o=y2o,
##                d2h=10,d2d=14,IFR=propd,HFR=proph)
## ## varying outrate version
## dorun <- function(x){
##   y[1:lkdnpt] <- 1
##   y[(1+lkdnpt):length(y)] <- exp(x[3])
##   md <- seirid2(I0=exp(x[1]),R0=exp(x[2]),Rinit = Rinit/UR,y=y,tt=tz,
##                d2h=delays$dC2H,d2d=delays$dC2D,y2o=y2o,
##                IFR=propd,HFR=exp(x[4]))
##   md$run(tz)
## }

seiridLL <- function(x,printbits=FALSE){
  out <- dorun(x)
  nd <- nrow(tgts)
  ecss <- out[1:nd,'incidence']
  disch <- out[1:nd,'disch']
  admns <- out[1:nd,'admns']
  dths <- out[1:nd,'deaths']
  LL1 <- sum(dpois(tgts$cases,UR*ecss,log=TRUE),na.rm=TRUE)  #cases
  LL2 <- sum(dpois(tgts$deaths,dths,log=TRUE),na.rm=TRUE)  #deaths
  LL3 <- sum(dpois(tgts$newhosp,admns,log=TRUE),na.rm=TRUE)  #admission
  LL4 <- 0#sum(dpois(tgts$disch,disch,log=TRUE),na.rm=TRUE)  #discharges
  LL <- LL1 + LL2 + LL3 + LL4
  if(printbits) print(c(LL1,LL2,LL3,LL4))
  LL
}

## test
seiridLL(c(0,1,-0.1,-2),printbits = TRUE)                    #test

## compare sheffield
UR <- 0.075*1
resi <- optim(par=c(0,1,-0.1,-2)+runif(4)/100,fn=seiridLL,
              control = list(fnscale=-1),hessian=TRUE) #ML
(Rzero <- exp(resi$par[2]))
(Effect <- 100*(1-exp(resi$par[3])))
(Rnet <- exp(sum(resi$par[2:3])))
(HFR <- 1e2*exp(sum(resi$par[4])))

cat(Rzero,file=here::here('data/Rzero.txt'))
cat(Effect,d=1,file=here::here('data/Effect.txt'))
cat(Rnet,file=here::here('data/Rnet.txt'))
cat(HFR,file=here::here('data/HFR.txt'))
cat(round(1e2*proph),file=here::here('data/HFRlit.txt'))


## compare sheffield
UR <- 0.075/2
resi2 <- optim(par=c(0,1,-0.1,-2)+runif(4)/100,fn=seiridLL,
              control = list(fnscale=-1),hessian=TRUE) #ML
(Rzero <- exp(resi$par[2]))
(Effect <- 100*(1-exp(resi$par[3])))
(Rnet <- exp(sum(resi$par[2:3])))
(HFR <- 1e2*exp(sum(resi$par[4])))
UR2 <- UR
UR <- UR*2

xx <- resi$par
seiridLL(xx)

## make dates for the outputs
dtz <- rep(dmy(lkdsig),length(tz))
shift <- seq(from=-24,by=1,length.out=length(tz))
which(shift==0)
dtz <- dtz + days(shift)
tmap <- data.table(t=tz,date=dtz)

## MLE for basecase as test
outi <- dorun(resi$par)
LA <- as.data.table(outi)
LA[,confirmed:='7.5%']
LA[,id:=0.0]
LA[,dys:=dtz]
ou1 <- LA
ou1[,trueincidence:=incidence]
ou1[,notes:=incidence*UR]
ou1[,ccadm:=propcc*trueincidence]        #to change
ou1[,beds:=hosp]
ou1[,hosp:=admns]
US1 <- ou1[,.(date=dys,confirmed,id,
            cases=notes,
            truecases=trueincidence,
            hosp,ccadm,deaths,
            dys)]
## reformat
SU1 <- melt(US1[,.(date=(date),confirmed,id,
                  cases,truecases,hosp,ccadm,deaths
                  )],id.vars = c('date','confirmed','id'))
names(SU1)[4] <- 'quantity'
SUM1 <- SU1[,.(mid=value,hi=value,lo=value),by=.(date,quantity,confirmed)]


## checks and uncertainty
## this is rather an under-estimate
## uncertainty over parameter estimates
Sig <- solve(-resi$hessian)
Sig2 <- solve(-resi2$hessian)
PMZ <- mvrnorm(200,resi$par,Sigma=Sig)
PMZ2 <- mvrnorm(200,resi2$par,Sigma=Sig2)

## CIs
## basecase
Rzero <- exp(PMZ[,2])
(Rzero <- c(mean=mean(Rzero),quantile(Rzero,.025),quantile(Rzero,.975)))
Rzero <- round(Rzero,digits=1)
Rzero <- paste0((Rzero[1])," (",Rzero[2],' to ',Rzero[3],")")
cat(Rzero,file=here::here('data/Rzero.txt'))
Effect <- 100*(1-exp(PMZ[,3]))
(Effect <- c(mean=mean(Effect),quantile(Effect,.025),quantile(Effect,.975)))
Effect <- round(Effect,digits=0)
Effect <- paste0((Effect[1])," (",Effect[2],' to ',Effect[3],")")
cat(Effect,file=here::here('data/Effect.txt'))
Rnet <- exp(rowSums(PMZ[,2:3]))
(Rnet <- c(mean=mean(Rnet),quantile(Rnet,.025),quantile(Rnet,.975)))
Rnet <- round(Rnet,digits=1)
Rnet <- paste0((Rnet[1])," (",Rnet[2],' to ',Rnet[3],")")
cat(Rnet,file=here::here('data/Rnet.txt'))
HFR <- exp(PMZ[,4])
(HFR <- 100*c(mean=mean(HFR),quantile(HFR,.025),quantile(HFR,.975)))
HFR <- round(HFR,digits=1)
HFR <- paste0((HFR[1])," (",HFR[2],' to ',HFR[3],")")
cat(HFR,file=here::here('data/HFR.txt'))

## SA
Rzero <- exp(PMZ2[,2])
(Rzero <- c(mean=mean(Rzero),quantile(Rzero,.025),quantile(Rzero,.975)))
Rzero <- round(Rzero,digits=1)
Rzero <- paste0((Rzero[1])," (",Rzero[2],' to ',Rzero[3],")")
cat(Rzero,file=here::here('data/Rzero2.txt'))
Effect <- 100*(1-exp(PMZ2[,3]))
(Effect <- c(mean=mean(Effect),quantile(Effect,.025),quantile(Effect,.975)))
Effect <- round(Effect,digits=0)
Effect <- paste0((Effect[1])," (",Effect[2],' to ',Effect[3],")")
cat(Effect,file=here::here('data/Effect2.txt'))
Rnet <- exp(rowSums(PMZ2[,2:3]))
(Rnet <- c(mean=mean(Rnet),quantile(Rnet,.025),quantile(Rnet,.975)))
Rnet <- round(Rnet,digits=1)
Rnet <- paste0((Rnet[1])," (",Rnet[2],' to ',Rnet[3],")")
cat(Rnet,file=here::here('data/Rnet2.txt'))
HFR <- exp(PMZ2[,4])
(HFR <- 100*c(mean=mean(HFR),quantile(HFR,.025),quantile(HFR,.975)))
HFR <- round(HFR,digits=1)
HFR <- paste0((HFR[1])," (",HFR[2],' to ',HFR[3],")")
cat(HFR,file=here::here('data/HFR2.txt'))


LU2 <- LU <- list()
for(i in 1:nrow(PMZ)){
  LU[[i]] <- as.data.table(dorun(PMZ[i,]))
  LU[[i]][,id:=i]
  UR <- UR/2
  LU2[[i]] <- as.data.table(dorun(PMZ2[i,]))
  LU2[[i]][,id:=i]
  UR <- UR*2
}

LU <- rbindlist(LU);  LU2 <- rbindlist(LU2)
LU[,confirmed:='7.5%'];  LU2[,confirmed:='3.8%']
LU <- rbind(LU,LU2)

ou <- merge(tmap,LU,by.x='t',by.y='t')

## version with multiple runs
ou[,trueincidence:=incidence]
ou[,notes:=incidence*UR]
ou[,ccadm:=propcc*trueincidence]        #to change
ou[,beds:=hosp]
ou[,hosp:=admns]
ou[confirmed=="3.8%",notes:=incidence*UR2]

US <- ou[,.(date,confirmed,id,
            cases=notes,
            truecases=trueincidence,
            hosp,ccadm,deaths)]
            ## dys)]

## reformat
SU <- melt(US[,.(date=(date),confirmed,id,
                  cases,truecases,hosp,ccadm,deaths
                  )],id.vars = c('date','confirmed','id'))
names(SU)[4] <- 'quantity'


SUM <- SU[,.(mid=mean(value),#quantile(value,0.5),
             hi=quantile(value,0.975),
             lo=quantile(value,0.025)),
          by=.(date,quantity,confirmed)]

## graphing targets
inc <- tgts[,.(date,cases,deaths,hosp=newhosp)]
inc <- melt(inc,id = 'date')
inc[,c('mid','confirmed','quantity','hi','lo'):=list(value,NA,variable,NA,NA)]

GPi <- ggplot(inc,aes(date,mid,col=quantity,label=mid)) + geom_point() +
  geom_line() + scale_y_sqrt() +
  geom_text_repel(show.legend = FALSE)+
  xlab('Date') + ylab('Daily incidence of quantity (square root scale)') + 
  theme(legend.position = c(0.15, 0.75),legend.direction='vertical')

pnm <- glue('plots') + '/2incData_' + td + '.pdf'       #
ggsave(GPi,filename = pnm,w=7,h=5)
ggsave(GPi,filename = here::here('figs/IncData.pdf'),w=7,h=5)


## incidence plot
GP3ud <- ggplot(SUM[date<=dmy('15/07/2020') & !quantity %in% c('ccadm','truecases')],
              aes(date,mid,col=quantity,lty=confirmed)) +
  geom_point(data=inc) + 
  geom_line() +
  geom_ribbon(aes(ymin=lo,ymax=hi,fill=confirmed),alpha=.2,col=NA)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  xlab('Date') + ylab('Daily incidence of quantity') +
  scale_y_continuous(label=comma) + 
  scale_color_manual(breaks=c('truecases','cases','hosp','ccadm','deaths'),
                     labels=c('true incidence','confirmed cases','hospital admissions',
                              'critical care admissions','deaths'),
                     values=cbbPalette[1:5])+
  scale_linetype_manual(breaks=c("7.5%","3.8%"),values=2:1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(col=guide_legend(ncol=2))+ guides(fill=FALSE)+
  theme(legend.position = c(0.75, 0.15),legend.direction='vertical') +
  geom_vline(xintercept = ymd(pnts$date[lkdnpt]),col=2,lty=2)+
  ggtitle('Projected daily numbers for Sheffield\nReal scale with 95% CIs') +
  facet_wrap(scales='free_y',~quantity,ncol=2)
GP3ud

GP3udL <- GP3ud + scale_y_log10(label=comma) +
  ylab('Daily incidence of quantity (log scale)') +
  ggtitle('Projected daily numbers for Sheffield\nReal scale with 95% CIs')
GP3udL


pnm <- glue('plots') + '/2incU_' + td + '.pdf'       #
ggsave(GP3ud,filename = pnm)
pnm <- glue('plots') + '/2LincU_' + td + '.pdf'       #
ggsave(GP3udL,filename = pnm)

ggsave(GP3ud,filename = here::here('figs/IncReal.pdf'))
ggsave(GP3udL,filename = here::here('figs/IncLog.pdf'))


## --- prevs
US <- merge(tmap,LU[,.(t,id,hosp,confirmed)],by.x='t',by.y='t')
US2 <- merge(tmap,LU2[,.(t,id,hosp,confirmed)],by.x='t',by.y='t') #jj

US <- rbind(US,US2)
US[,HDU:=hosp*bedprops[variable=='Number of confirmed COVID 19 HDU patients',prop]]
US[,ITU:=hosp*bedprops[variable=='Number of confirmed COVID 19 ITU patients',prop]]
US[,IDU:=hosp*bedprops[variable==vrs[3],prop]]
US[,other:=hosp*bedprops[variable=='Number of confirmed COVID 19 any other beds',prop]]

US[,O2:=hosp*o2props[variable=='Number of confirmed COVID 19 on O2',prop]]
US[,NIV:=hosp*
      o2props[variable=='Number of confirmed COVID 19 on Non Invasive Ventilation',prop]]
US[,Mechanical:=hosp*
      o2props[variable=='Number of confirmed COVID 19 on Mechanical Ventilation',prop]]

## US[,confirmed:=NA]
## nb 33% national vs 12%
names(US)
## US[,c('dys','cases','Population','dta','ddys'):=NULL]

## reformat
PU <- melt(US,id.vars = c('date','confirmed','id'))
names(PU)[4] <- 'quantity'

PUM <- PU[,.(mid=mean(value),#quantile(value,0.5),
             hi=quantile(value,0.975),
             lo=quantile(value,0.025)),
          by=.(date,quantity,confirmed)]

## cc & hosp etc
IDUP <- DM[variable==vrs[3],.(mid=value),by=date] #IDU
IDUP[,quantity:='IDU']
HDUP <- DM[variable==vrs[1],.(mid=value),by=date] #HDU
HDUP[,quantity:='HDU']
MP <- DM[variable==vrs[7],.(mid=value),by=date] #MV
MP[,quantity:='Mechanical']
NP <- DM[variable==vrs[6],.(mid=value),by=date] #NIV
NP[,quantity:='NIV']
OP <- DM[variable==vrs[5],.(mid=value),by=date] #O2
OP[,quantity:='O2']
RP <- DM[variable==vrs[4],.(mid=value),by=date] #other beds
RP[,quantity:='other']
hosp <- DM[variable %in% vrs[1:4],.(mid=sum(value)),by=date]
hosp[,quantity:='hosp']
ccadm <- DM[variable == vrs[2],.(mid=(value)),by=date] #NOTE CC = ITU assumed
ccadm[,quantity:='ITU']
prev <- rbindlist(list(IDUP,HDUP,MP,NP,OP,RP,hosp,ccadm))
prev[,confirmed:=NA]

PUM <- PUM[quantity!='t']



## prevalence plot
GP4ud <- ggplot(PUM[date<=dmy('15/07/2020') ],
              aes(date,mid,col=quantity,lty=confirmed)) +
  geom_line() +
  geom_point(data=prev) +
  geom_ribbon(aes(ymin=lo,ymax=hi,fill=confirmed),alpha=.2,col=NA)+  
  xlab('Date') + ylab('Prevalence of quantity') +
  scale_y_continuous(label=comma) +
  scale_color_manual(values=cbbPalette[1:8])+
  scale_linetype_manual(breaks=c("7.5%","3.8%"),values=2:1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(col=guide_legend(ncol=4))+ guides(fill=FALSE)+
  guides(linetype=guide_legend(ncol=2)) + 
  theme(legend.position = c(0.85, 0.1),legend.direction='vertical') +
  ggtitle('Projected prevalent  numbers for Sheffield\nReal scale with 95% CIs') +
  expand_limits(x=PUM[,min(date)]) +
  facet_wrap(scales='free_y',~quantity)
GP4ud

GP4udL <- GP4ud + scale_y_log10(label=comma) +
  ylab('Prevalence of quantity (log scale)')  +
  ggtitle('Projected prevalent  numbers for Sheffield\nLog scale with 95% CIs') 
GP4udL


pnm <- glue('plots') + '/2PrevU_' + td + '.pdf'       #
ggsave(GP4ud,filename = pnm)
pnm <- glue('plots') + '/2LPrevU_' + td + '.pdf'       #
ggsave(GP4udL,filename = pnm)


ggsave(GP4ud,filename = here::here('figs/PrevReal.pdf'))
ggsave(GP4udL,filename = here::here('figs/PrevLog.pdf'))


## save most recent dates
cat(as.character(max(pnts$date)),file=here::here('data/maxdate_cases.txt'))
cat(as.character(max(DM$date)),file=here::here('data/maxdate_sitrep.txt'))


## save output also
fn <- glue(here::here('data'))+ '/PUM_' + td + '.Rdata'
save(PUM,file=fn)
fn <- glue(here::here('data'))+ '/SUM_' + td + '.Rdata'
save(SUM,file=fn)
SUM[quantity=='hosp',quantity:='admissions']
fwrite(PUM[date<=dmy('15/07/2020') ],file=here::here('data/PrevData.csv'))
fwrite(SUM[date<=dmy('15/07/2020') & !quantity %in% c('ccadm','truecases')],
       file=here::here('data/IncData.csv'))
nmz <- c('ages','demo','sympto','symptohosp','hospcc','IFR')
fwrite(format(parms[,..nmz],digits=1,scientific=FALSE),
       file=here::here('data/parms.csv'))

