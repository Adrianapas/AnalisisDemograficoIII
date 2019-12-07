##################################
####  Proyección Nuevo León  #####
##################################
####################################################
#### Maestría en Demografía                    #####
#### El Colegio de México                      #####
#### Autores:  Adriana Perez-Arciniega Soberon #####
####           Isaac Díaz Martínez             #####
#### Fecha: 06/12/2019                         #####
####################################################

### Librerias

rm(list=ls())
gc()
require(forecast)
require(ggplot2)
require(mvtnorm)
require(simPop)

### Mortalidad

mxNL=read.csv("mx.csv",header=T)
mxNL<-mxNL[,-c(1,2)]

### Migración Internacional

InmxNL=read.csv("Imx5.csv",header=T)
EmigxNL=read.csv("Emx5.csv",header=T)
NxNL=read.csv("Nx.csv",header=T)

InmxNL.F5<-data.frame(InmxNL[InmxNL$Sexo=="Mujer",-1])
InmxNL.H5<-data.frame(InmxNL[InmxNL$Sexo=="Hombre",-1])
EmigxNL.F5<-data.frame(EmigxNL[EmigxNL$Sexo=="Mujer",-1])
EmigxNL.H5<-data.frame(EmigxNL[EmigxNL$Sexo=="Hombre",-1])

InmxNL.F=matrix(0,81,46)
InmxNL.H=matrix(0,81,46)
EmigxNL.F=matrix(0,81,46)
EmigxNL.H=matrix(0,81,46)

for(i in 0:45){
  InmxNL.F[,i+1]=as.matrix(sprague(InmxNL.F5[,i+2]))
  InmxNL.H[,i+1]=as.matrix(sprague(InmxNL.H5[,i+2]))
  EmigxNL.F[,i+1]=as.matrix(sprague(EmigxNL.F5[,i+2]))
  EmigxNL.H[,i+1]=as.matrix(sprague(EmigxNL.H5[,i+2]))
}

maNxM=as.matrix(NxNL[NxNL$Sexo=="Hombres",-c(1:2)])
maNxF=as.matrix(NxNL[NxNL$Sexo=="Mujeres",-c(1:2)])

ixt.MNL<-InmxNL.H/maNxM
ixt.FNL<-InmxNL.F/maNxF
ext.MNL<-EmigxNL.H/maNxM
ext.FNL<-EmigxNL.F/maNxF

### Fecundidad

fx5<- read.csv("Fecundidad_NL.CSV",header=TRUE)
fx5q<- read.csv("fxVer.CSV",header=TRUE)

asfr=function(fx5,year){
  
  fx0=fx5[fx5$Year==year,"fx"]
  
  Fx=5*cumsum(fx0)
  TGF=Fx[7]
  FxF=Fx/TGF
  
  x5=seq(17.5,47.5,5)
  Yx= log(-log(FxF))
  
  Yx.lm=lm(Yx[-7] ~ x5[-7])
  a=Yx.lm$coefficients[1]
  b=Yx.lm$coefficients[2]
  
  A=exp(-exp(a))
  B=exp(b)
  x1=c(15:50)
  
  (Fx.estim=TGF*A^(B^(x1)))
  
  fx=Fx.estim[2:36]-Fx.estim[1:35]
  fx=c(Fx.estim[1],fx)
  
  return(fx)
}

fx1=data.frame(matrix(0,36,46))
row.names(fx1)=c(15:50)
names(fx1)=c(1970:2015)

for(i in 1:46){
  
  fx1[,i]=asfr(fx5q,1969+i)
  
}

fx1=fx1[1:35,19:45]


### Elementos

edades <- dim(mxNL)[1] #dame la dimension y quedate con el numero de renglones = 220
tiempo.mort <- dim(mxNL)[2] #quedo con el número de columnas
añoini.mort <- 1970
añobase <- 2015
horizonte <- 25
añofin <- añobase+horizonte
tiempo.tot <- tiempo.mort+horizonte
edades.fec <-dim(fx1)[1]
añoini.fec <- 1995
tiempo.fec <-dim(fx1)[2]
edades.migNL <- dim(ixt.FNL)[1]
tiempo.migNL <- dim(ixt.FNL)[2]
añoini.migNL <- 1995

### Inicia la estimacion con el método de Lee Carter (1992)


lc.svd <- function(m,edades,tiempo1, tiempo2,ln){
  if(ln==TRUE){
    lm <- log(m)
  }else{
    lm<-
      m
  }
  
  ax <- rowMeans(lm[,tiempo1:tiempo2])
  lm_a <- lm -ax
  d <- matrix(0,nr = min(edades,tiempo2),
              nc = min(edades,tiempo2))
  
  diag(d) <- svd(lm_a)$d
  
  kt <-(d%*%t(-svd(lm_a)$v))
  bx <- -svd(lm_a)$u
  
  lc.svd=list(ax = ax, bx = bx, kt = kt, D = d)
  
  
}
tabmort <- function(m,edades,sex){
  
  mx <- m
  
  nax <- matrix(0.5,dim(mx)[1],dim(mx)[2])
  ## 1 MUJERES 2 HOMBRES
  if(sex==1){
    for(i in 1:dim(mx)[2]){
      if(mx[1,i]<0.01724){
        nax[1,i] <- 0.14903-2.05527*mx[1,i]
      }else if(mx[1,i]>=0.01724 & mx[1,i]<0.06891){
        nax[1,i] <- 0.04667+3.88089*mx[1,i]
      }else{nax[1,i] <- 0.31411}
    }
  }else{
    for(i in 1:dim(mx)[2]){
      if(mx[1,i]<0.023){
        nax[1,i] <- 0.14929-1.99545*mx[1,i]
      }else if(mx[1,i]>=0.023 & mx[1,i]<0.08307){
        nax[1,i] <- 0.02832+3.26021*mx[1,i]
      }else{nax[1,i] <- 0.29915}
    }
  }
  
  
  nax[edades,] <- 1/mx[edades,]
  
  qx<-matrix(1,dim(mx)[1],dim(mx)[2])
  
  for(i in 1:(dim(mx)[1])){
    qx[i,]<-mx[i,]/(1+(1-nax[i,])*mx[i,])
  }
  
  px <- 1-qx
  
  lx<-matrix(1,dim(mx)[1],dim(mx)[2])
  
  for(i in 2:dim(mx)[1]){
    lx[i,] <- lx[i-1,]*px[i-1,]
  }
  
  dx <- matrix(0,dim(mx)[1],dim(mx)[2])
  dx[dim(mx)[1],] <- lx[dim(mx)[1],]
  for(i in 1:(dim(mx)[1]-1)){
    dx[i,]<-lx[i,]-lx[i+1,]
  }
  
  
  Lx<-matrix(0,dim(mx)[1],dim(mx)[2])
  Lx[1,] <- dx[1,]/mx[1,]
  Lx[edades,] <- dx[edades,]/mx[edades,]
  for(i in 2:(edades-1)){
    Lx[i,]<-(lx[i,]+lx[i+1,])/2
  }
  
  Tx<-matrix(0,dim(mx)[1],dim(mx)[2])
  Tx[edades,]<-Lx[edades,]
  for(i in (edades-1):1){
    Tx[i,]<-Lx[i,]+Tx[i+1,]
  }
  
  ex <- Tx/lx
  
  Sx<-matrix(0,(dim(mx)[1]+1),dim(mx)[2])
  Sx[1,]<-Lx[1,]/lx[1,]
  Sx[(edades+1),] <- Tx[edades,]/Tx[(edades-1),]
  for(i in 2:edades){
    Sx[i,]<-Lx[i,]/Lx[i-1,]
  }
  
  tabmort <- list(Edad=c(0:(edades-1)),mx=mx, nax=nax, qx=qx, 
                  px=px, lx=lx, dx=dx, Lx=Lx, Tx=Tx, ex=ex, Sx=Sx)
}

lc.mort <-lc.svd(mxNL, edades, tiempo1 = 41, 
                 tiempo2 = tiempo.mort, 
                 ln=TRUE)

lc.fec <-lc.svd(fx1, edades=edades.fec, 
                tiempo1 = 10, 
                tiempo2 = tiempo.fec,
                ln=TRUE)

lc.inmFNL <-lc.svd(ixt.FNL, edades=edades.migNL, 
                   tiempo1 = 1, 
                   tiempo2 = tiempo.migNL,
                   ln=TRUE)

lc.inmMNL <-lc.svd(ixt.MNL, edades=edades.migNL, 
                   tiempo1 = 1, 
                   tiempo2 = tiempo.migNL,
                   ln=TRUE)

lc.emigFNL <-lc.svd(ext.FNL, edades=edades.migNL, 
                  tiempo1 = 1, 
                  tiempo2 = tiempo.migNL,
                  ln=TRUE)

lc.emigMNL <-lc.svd(ext.MNL, edades=edades.migNL, 
                  tiempo1 = 1, 
                  tiempo2 = tiempo.migNL,
                  ln=TRUE)

kt1.fit <-auto.arima(lc.mort$kt[1,], trace=T, d=1)

ft1.fit <-auto.arima(lc.fec$kt[1,], trace=T, d=1)
                     
it1F.fit <- auto.arima(lc.inmFNL$kt[1,], trace=T, allowdrift = F)

it1M.fit <- auto.arima(lc.inmMNL$kt[1,], trace=T, allowdrift = F)

et1F.fit <- auto.arima(lc.emigFNL$kt[1,], trace=T, allowdrift = F,d=1)

et1M.fit <- auto.arima(lc.emigMNL$kt[1,], trace=T, allowdrift = F,d=1)


### Proyeccion

h <- 25

kt.for <- forecast(kt1.fit, h = h, c(95))

ft.for <- forecast(ft1.fit, h = h, c(95))

itF.for <- forecast(it1F.fit, h = h, c(95))
itM.for <- forecast(it1M.fit, h = h, c(95))

etF.for <- forecast(et1F.fit, h = h, c(95))
etM.for <- forecast(et1M.fit, h = h, c(95))


### Tasas centrales de mortalidad
mx.for <- exp(lc.mort$ax + lc.mort$bx[,1]%*%t(kt.for$mean))

### Funciones de supervivencia 
SxF.for <- tabmort(mx.for[111:220,], edades = 110, sex=1)$Sx
SxM.for <- tabmort(mx.for[1:110,], edades = 110, sex=2)$Sx

### Tasas especificas de fecundidad
fx.for<-exp(lc.fec$ax + lc.fec$bx[,1]%*%t(ft.for$mean))
TGF=colSums(fx.for)           

### Tasas de migración
ixF.For <- rbind(exp(lc.inmFNL$ax + lc.inmFNL$bx[,1]%*%t(itF.for$mean)),
                 matrix(0,0,25))

ixM.For <- rbind(exp(lc.inmMNL$ax + lc.inmMNL$bx[,1]%*%t(itM.for$mean)),
                 matrix(0,0,25))

exF.For <- rbind(exp(lc.emigFNL$ax + lc.emigFNL$bx[,1]%*%t(etF.for$mean)),
                 matrix(0,0,25))

exM.For <- rbind(exp(lc.emigMNL$ax + lc.emigMNL$bx[,1]%*%t(etM.for$mean)),
                 matrix(0,0,25))

### Población

Px=read.csv("Px.csv",header=T)

PxFNL <- Px[Px$Sexo=="Mujeres", -c(1,2)]
PxMNL <- Px[Px$Sexo=="Hombres", -c(1,2)]

NxFNL <- Px[NxNL$Sexo=="Mujeres", -c(1,2)]
NxMNL <- Px[NxNL$Sexo=="Hombres", -c(1,2)]

PxF.for <- matrix(0,81,26)
PxM.for <- matrix(0,81,26)

NxF.for <- matrix(0,81,26)
NxM.for <- matrix(0,81,26)

PxF.for[,1] <- PxFNL[,"X2016"]
PxM.for[,1] <- PxMNL[,"X2016"]

NxF.for[,1] <- NxFNL[,"X2015"]
NxM.for[,1] <- NxMNL[,"X2015"]

Bx <- matrix(0,35,26)
BF <- vector(length=26)
BM <- vector(length=26)

# MUJERES 
for(i in 2:26){
  
  #Con pob a mitad de año 1 a 108
  PxF.for[2:80,i] <- (PxF.for[1:79,i-1] + 
                         0.5*NxF.for[1:79,i-1]*ixF.For[1:79,i-1])* SxF.for[1:79,i-1]+
    NxF.for[2:80,i-1]*0.5*ixF.For[2:80,i-1]-
    NxF.for[1:79,i-1]*exF.For[1:79,i-1]
  
  #ultimo grupo 109 y mas
  PxF.for[81,i]<-(PxF.for[80,i-1] + 
                     0.5*NxF.for[80,i-1]*ixF.For[80,i-1])*SxF.for[80,i-1] -
    NxF.for[80, i-1]*exF.For[80,i-1] +
    (PxF.for[81,i-1] + 
       NxF.for[81,i-1]*0.5*ixF.For[81,i-1])*SxF.for[81,i-1]+
    NxF.for[81,i-1]*0.5*ixF.For[81,i-1]-
    NxF.for[81,i-1]*0.5*exF.For[81,i-1]
  
  
  #Nacimientos 
  Bx[,i-1] <- fx.for[,i-1]*(PxF.for[16:50,i-1]+
                              0.5*NxF.for[16:50,i-1]*ixF.For[16:50,i-1] +
                              PxF.for[16:50,i]) * 0.5
  
  BF[i-1] <- (1/2.05)*sum(Bx[,i-1])
  
  
  #primer grupo de edad
  
  PxF.for[1,i] <- BF[1]*SxF.for[1,i-1] + 
    NxF.for[1,i-1]*0.5*ixF.For[1,i-1]+
    NxF.for[1,i-1]*exF.For[1,i-1]
  
  
  #POb mitad de año
  
  NxF.for[,i] <- 0.5*(PxF.for[,i-1] + PxF.for[,i])
}

# HOMBRES
for(i in 2:26){
  
  #Con pob a mitad de año 1 a 108
  PxM.for[2:80,i] <- (PxM.for[1:79,i-1] + 
                         0.5*NxM.for[1:79,i-1]*ixM.For[1:79,i-1])* SxM.for[1:79,i-1]+
    NxM.for[2:80,i-1]*0.5*ixM.For[2:80,i-1]-
    NxM.for[1:79,i-1]*exM.For[1:79,i-1]

  #ultimo grupo 109 y mas
  PxM.for[81,i]<-(PxM.for[80,i-1] + 
                     0.5*NxM.for[80,i-1]*ixM.For[80,i-1])*SxM.for[80,i-1] -
    NxM.for[80, i-1]*exM.For[80,i-1] +
    (PxM.for[81,i-1] + 
       NxM.for[81,i-1]*0.5*ixM.For[81,i-1])*SxM.for[81,i-1]+
    NxM.for[81,i-1]*0.5*ixM.For[81,i-1]-
    NxM.for[81,i-1]*0.5*exM.For[81,i-1]
  
  
  #Nacimientos 
  
  BM[i-1] <- (1.05/2.05)*sum(Bx[,i-1])
  
  
  #primer grupo de edad
  
  PxM.for[1,i] <- BM[1]*SxM.for[1,i-1] + 
    NxM.for[1,i-1]*0.5*ixM.For[1,i-1]+
    NxM.for[1,i-1]*exM.For[1,i-1]
  
  
  #POb mitad de año
  
  NxM.for[,i] <- 0.5*(PxM.for[,i-1] + PxM.for[,i])
}

### Población total proyectada
POB <- colSums(PxF.for) + colSums(PxM.for)          
POB

### Graficas

#Mortalidad

matplot(cbind(log(mx.for[1:110,1]),log(mx.for[111:220,1]),log(mx.for[1:110,25]),log(mx.for[111:220,25])),type="l",xlab="Edades desagrupadas",ylab="Log(Tasas centrales de mortalidad)",main="Proyección de las tasas centrales de mortalidad",sub="Estado: Nuevo León",lwd=3,col = c("blue","red","lightgoldenrod4","pink" ))
legend(0,-0.5,c("Hombres 2015","Mujeres 2015","Hombres 2040","Mujeres 2040"), 
       lty = c(1,2),col = c("blue","red","lightgoldenrod4","pink"),cex = 1)
grid(col="orange")

#Migración

ixF.For1=ixF.For*NxF.for[,-26]
exF.For1=exF.For*NxF.for[,-26]
SNMF=ixF.For1-exF.For1

ixM.For1=ixM.For*NxF.for[,-26]
exM.For1=exM.For*NxF.for[,-26]
SNMM=ixM.For1-exM.For1

SNMT=SNMM+SNMF
sNMtota=colSums(SNMT)

inmig.tot=ixF.For1+ixM.For1
inmig.total=colSums(inmig.tot)

emig.tot=exF.For1+exM.For1
emig.total=colSums(emig.tot)

años=seq(2016,2040,1)
matplot(años,cbind(sNMtota,inmig.total,emig.total),col = c("blue","lightblue","red"),type="l",ylab="Stock de migrantes",xlab="Años proyectados",main="Saldo neto migratorio internacional proyectado",sub="Estado: Nuevo León (2015:2040)",lwd=3)
legend(2016,2500,c("Saldo neto migratorio internacional","Inmigración internacional","Emigración internacional"),
       lty = c(1,2),col = c("blue","lightblue","red"),cex = 0.8)
grid(col="orange")

#Fecundidad

fx.forup<-exp(lc.fec$ax + lc.fec$bx[,1]%*%t(ft.for$upper))
TGFup=colSums(fx.forup)           
fx.forlow<-exp(lc.fec$ax + lc.fec$bx[,1]%*%t(ft.for$lower))
TGFlowe=colSums(fx.forlow)           

par(mfrow=c(1,2))
años=seq(2016,2040,1)
matplot(años,cbind(TGF,TGFup,TGFlowe),col = c("blue","lightblue","red"),type="l",ylab="Promedio de hijos por mujer",xlab="Años proyectados",main="Proyección de tasa global de fecundidad (TGF)",sub="Estado: Nuevo León (2015:2040)",lwd=3)
legend(2016,2.79,c("Escenario Medio","Escenario Alto","Escenario Bajo"),
       lty = c(1,2),col = c("blue","lightblue","red"),cex = 0.8)
grid(col="orange")

Edades=seq(15,49,1)
matplot(Edades,fx.for,type="l",ylab="fx",xlab="Edades desagrupadas",main="Proyección de las tasas específicas de fecundidad",sub="Estado: Nuevo León (2015:2040)",lwd=2)
grid(col="orange")

#Población Total

POBF=colSums(PxF.for)
POBM=colSums(PxM.for)
POB = POBF+POBM         

años=seq(2015,2040,1)
par(mfrow=c(1,2))
matplot(años,cbind(POB,POBF,POBM),col = c("blue","lightblue","red"),type="l",ylab="Poblacion",xlab="Años proyectados",main="Población total proyectada",sub="Estado: Nuevo León (2015:2040)",lwd=3)
legend(2015.5,5400000,c("Población total","Población femenina","Población masculina"),
       lty = c(1,2),col = c("blue","lightblue","red"),cex = 0.8)
grid(col="orange")

POB.edad <- PxF.for+ PxM.for          
matplot(POB.edad[-81,],type="l",ylab="Poblacion",xlab="Edades",main="Estructura por edades de ambos sexos proyectada",sub="Estado: Nuevo León (2015:2040)",lwd=1)
grid(col="orange")
