
# basic function to plot a transiogram
tr_2d_vert<-function(datai,dataj,size,ncell,dz){
  cov<-numeric(size)
  for (l in 1:size){
    npair = nrow(datai)-l
    t =0
    for ( pair in 1:npair){
      for ( y in 1:ncol(datai)){
        if(is.na(datai[pair,y])== FALSE & is.na(dataj[pair+l,y]) == FALSE){
          cov[l] = cov[l]+ datai[pair,y]*dataj[pair+l-1,y]
          t = t+1
        }
      }
    }
    cov[l]=cov[l]/t
  }
  cov = cov /mean(datai,na.rm = TRUE)# we are measuring a correlation
  dist = (0:(size-1))/(ncell*dz)
  plot(dist,cov,type = "p",col="black",pch=21,bg="grey",cex=2,ylim=c(0.,1),ylab = "", xlab = "", axes = FALSE)
  axis(1,at = seq(0,3,0.5),cex.axis=1.5)
  axis(2,cex.axis=1.5)
  #mtext("Distance", side = 1, line = 2.5, cex = 1.)
  #mtext("Transiogram", side = 2, line = 2.5, cex = 1.,las = 0)
  return(cov)
}

embedded3d<-function(datai,dataj){
  # we assume there is no NA 
  ndim=dim(datai)
  embed=0
  Tii=0
  n=ndim[3]-1
  prop=mean(datai[,,-182]) # the last line is not counted
  for ( z in 1:n){
    Tii=Tii+mean(datai[,,z]*datai[,,z+1])/n
    embed=embed+mean(datai[,,z]*dataj[,,z+1])/n
  }
  Tii=(1-Tii/prop)
  embed=embed/prop
  return(embed/Tii) 
}

embedded2d<-function(datai,dataj){
  
  embed=0
  Tii=0
  n=nrow(datai)-1
  prop=mean(datai[-(n+1),],na.rm=TRUE) # the last line is not counted
  count=0
  for ( z in 1:n){
    for (x in 1:ncol(dataj)){
      if(is.na(datai[[z,x]])==FALSE & is.na(datai[[z+1,x]])==FALSE){
        Tii=Tii+datai[[z,x]]*datai[[z+1,x]]
        embed=embed+datai[[z,x]]*dataj[[z+1,x]]
        count=count+1
      }
      
    }
  }
  Tii=(1-Tii/(prop*count))
  embed=embed/(prop*count)
  return(embed/Tii) 
}

sim2d<-function(rv,rh,freq,N,dz,nz,dx,nx){
  phase <- runif(N,0,2*pi)
  bin <- sample(c(-1,1),N,replace=T)
  w1 <- rnorm(N, mean = freq , sd = sqrt(2)/rv)
  w2 <- rnorm(N, mean = 0 , sd = sqrt(2)/rh)
  sim<-array(0,c(nx,nz))
  for(i in 1:N){
    for ( z in 1:nz){
      for(x in 1:nx){
        sim[x,z] = sim[x,z] + cos(bin[i]*w1[i]*(z*dz)+w2[i]*x*dx+phase[i])
      }
    }
  }
  sim = sim *sqrt(2/N)
  
  return(sim)
}

sim_specgauss1D<-function(rv,freq,N,dz,nz){
  phase <- runif(N,0,2*pi)
  bin <- sample(c(-1,1),N,replace=T)
  w1 <- rnorm(N, mean = freq , sd = sqrt(2)/rv)
  sim<-numeric(nz)
  for(i in 1:N){
    for ( z in 1:nz){
        sim[z] = sim[z] + cos(bin[i]*w1[i]*(z*dz)+phase[i])
    }
  }
  sim = sim *sqrt(2/N)
  
  return(sim)
}

sim3d<-function(rv,rhx,rhy,freq,N,dz,nz,dx,nx,dy,ny){
  phase <- runif(N,0,2*pi)
  bin <- sample(c(-1,1),N,replace=T)
  w1 <- rnorm(N, mean = freq , sd = sqrt(2)/rv)
  w2 <- rnorm(N, mean = 0 , sd = sqrt(2)/rhx)
  w3 <- rnorm(N, mean = 0 , sd = sqrt(2)/rhy)
  sim<-array(0,c(nx,ny,nz))
  for(i in 1:N){
    for ( z in 1:nz){
      for(x in 1:nx){
        for(y in 1:ny){
        sim[x,y,z] = sim[x,y,z] + cos(bin[i]*w1[i]*(z*dz)+w2[i]*x*dx+w3[i]*y*dy+phase[i])
        }
      }
    }
  }
  sim = sim *sqrt(2/N)
  
  return(sim)
}

tr_3d_vert<-function(datai,dataj,size){
  cov<-numeric(size)
  ndim=dim(datai)
  for (l in 1:size){
    npair = ndim[3]-l
    for ( pair in 1:npair){
        cov[l] = cov[l]+ mean(datai[,,pair]*dataj[,,pair+l-1])/npair
    }
  }
  cov = cov /mean(datai)# we are measuring a correlation
  #dist = (0:(size-1))/(ncell*dz)
  #plot(dist,cov,type = "p",col="black",pch=21,bg="grey",cex=2,ylim=c(0.,1),ylab = "", xlab = "", axes = FALSE)
  #axis(1,at = seq(0,3,0.5),cex.axis=1.5)
  #axis(2,cex.axis=1.5)
  #mtext("Distance", side = 1, line = 2.5, cex = 1.)
  #mtext("Transiogram", side = 2, line = 2.5, cex = 1.,las = 0)
  return(cov)
}

sim2dexp<-function(rv,rh,freq,N,dz,nz,dx,nx){
  phase <- runif(N,0,2*pi)
  bin <- sample(c(-1,1),N,replace=T)
  w1 <- rcauchy(N, location = freq , scale =1/rv)
  w2 <- rnorm(N, mean = 0 , sd = sqrt(2)/rh)
  sim<-array(0,c(nx,nz))
  for(i in 1:N){
    for ( z in 1:nz){
      for ( x in 1:nx){
        sim[x,z] = sim[x,z] + cos(bin[i]*w1[i]*(z*dz)+w2[i]*x*dx+phase[i])
      }
    }
  }
  sim = sim *sqrt(2/N)
  
  return(sim)
}