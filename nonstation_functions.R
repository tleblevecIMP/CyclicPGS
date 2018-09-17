tr_2d_vert<-function(datai,dataj,size,ncell,dz){
  cov<-numeric(size)
  for (l in 1:size){
    npair = nrow(datai)-l
    t =0
    for ( pair in 1:npair){
      for ( y in 1:ncol(datai)){
        if(is.na(datai[pair,y])== FALSE & is.na(dataj[pair+l,y]) == FALSE){
          cov[l] = cov[l]+ datai[pair,y]*dataj[pair+l,y]
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
  mtext("Distance", side = 1, line = 2.5, cex = 1.5)
  #mtext("Transiogram", side = 2, line = 2.5, cex = 1.5,las = 0)
  return(cov)
}

embedded3d<-function(datai,dataj){
  # we assume there is no NA 
  ndim=dim(datai)
  embed=0
  Tii=0
  n=ndim[3]-1
  prop=mean(datai[,,-(n+1)]) # the last line is not counted
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
        if(is.na(datai[z,x])==FALSE & is.na(datai[z+1,x])==FALSE){
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


#ok, here we are coding the kriging of the proportions 

# we need to implement the surface kriging otherwise it is singular
prop_kriging<-function(grid,range,dx,dy,dz){

  data_loc<-which(!is.na(grid[,,1]),arr.ind=TRUE)
  estim_loc <-which(is.na(grid[,,1]),arr.ind=TRUE) #one surface
  #data_values=grid[!is.na(grid)]
  m = nrow(estim_loc)
  n=nrow(data_loc)
  #building the covariance matrix between the data points
  dist<-matrix(0,n,n)
  c_data<-matrix(0,n,n)
  for ( i in 1:(n-1)){
    for (j in (i+1):n ){
      distx=(data_loc[i,1]-data_loc[j,1])*dx
      disty=(data_loc[i,2]-data_loc[j,2])*dy
      dist[i,j]=dist[j,i]=sqrt((distx^2)+(disty^2))
    }
  }
  c_data<- exp(-(dist/range)^2) 
  invc<-chol2inv(chol(c_data)) #inverting the matrix

  # covariance vector between estimation and data of every surface
  # one should be enough to compute because all wells have same location
  c_estim_data<-matrix(0,m,n)
  for ( j in 1:m){
    distx=(data_loc[,1]-estim_loc[j,1])*dx
    disty=(data_loc[,2]-estim_loc[j,2])*dy
    dist=sqrt((distx^2)+(disty^2))
    c_estim_data[j,]<-exp(-(dist/range)^2)
  }
  
  # for every surface
  ndim<-dim(grid)
  for ( i in 1:ndim[3]){
    surf<-grid[,,i]
    data_values<-surf[!is.na(surf)]
    b <- t(data_values) %*% invc
    for ( j in 1:m){
      hx = estim_loc[j,1]
      hy = estim_loc[j,2]
      grid[hx,hy,i] = b %*% c_estim_data[j,]
      #lambdas=solve(c_data,c_estim_data[j,])
      #grid[hx,hy,i] = lambdas %*% data_values
    }
  }
  
  return(grid)
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

cov_data<-function(data_loc,dz,dhx,dhy,rhx,rhy,rv,b){ # coded for optimization
  n=nrow(data_loc)
  dist_hx <- dist_hy <- dist_v<- matrix(0,n,n)
  for ( i in 1:(n-1)){
    for (j in (i+1):n ){
      dist_v[i,j]=dist_v[j,i]=abs(data_loc[i,3]-data_loc[j,3])*dz
      dist_hx[i,j]=dist_hx[j,i]=(data_loc[i,1]-data_loc[j,1])*dhx # no abs because squared
      dist_hy[i,j]=dist_hy[j,i]=(data_loc[i,2]-data_loc[j,2])*dhy
    }
  }
  
  cov_h=exp(-sqrt( ((dist_hx/rhx)^2) + ((dist_hy/rhy)^2) )) #a verifier
  cov_v =exp(-(dist_v/rv)^2)*cos(b*dist_v)
  #separable model
  return(cov_h*cov_v)
}

negrnorm<-function(q,mean,sd){
  a=rnorm(1,mean,sd)
  while(a>q){
    a=rnorm(1,mean,sd)
  }
  return(a)
}

posrnorm<-function(q,mean,sd){
  a=rnorm(1,mean,sd)
  while(a<q){
    a=rnorm(1,mean,sd)
  }
  return(a)
}

truncnorm<-function(qmin,qmax,mean,sd){
  a=rnorm(1,mean,sd)
  while(a<qmin | a >qmax){
    a=rnorm(1,mean,sd)
  }
  return(a)
}

gibbs_sampling_map<-function(grid,qmap,dz,dhx,dhy,rv,rhx,rhy,b,it){
  # this algorithm is based on Emery paper (2013); algorithm 7
  
  data_values<-grid[!is.na(grid)] # getting the values
  q<-qmap[!is.na(qmap)]
  data_loc<-which(!is.na(grid),arr.ind=TRUE) # getting the position in the grid
  n=length(data_values)
  neg <-which(data_values==1)
  pos <- which(data_values==0)
  und <- which(data_values==Inf)
  # computation of the covariance matrix
  c_data<- cov_data(data_loc,dz,dhx,dhy,rhx,rhy,rv,b)
  inv <- which(c_data <0)
  #initialization uncorrelated gaussian simulation inside the interval
  undsim = rnorm(1,0,1)
  negsim= negrnorm(q[1],0,1) # we do this so that values do not start completely randomly
  possim= posrnorm(q[1],0,1)
  sim<-rep(0,n)
  for ( i in 1:n){
    if (data_values[i]==1){
      if(negsim>q[i]){
        negsim= negrnorm(q[i],0,1)
      }
      sim[i] = negsim
    }
    if (data_values[i]==0){ 
      if ( possim<q[i]){
        possim=posrnorm(q[i],0,1)
      }
      sim[i] = possim
      } 
    if (data_values[i]==Inf){sim[i]=undsim} # we might have a threshold here but it does not matter, no data
    #this can happen either if the facies is not constraining or if the well has been extended
  }
  
  # gibbs sampling
  for( t in 1:it){
    for ( i in 1:n){
      condi<-(((q-sim)/c_data[i,])+sim[i])
      maxi=Inf
      mini=-Inf
      for ( j in 1:n){
        # the covariance can be negative and so can change the sign of the inequality
        if (c_data[i,j]>0){
          if(data_values[j]==1){
            maxi = min(maxi,condi[j])
          }
          if(data_values[j]==0){
            mini=max(mini,condi[j])
            
          }
          # if data value = inf, no condition
        }
        if(c_data[i,j]<0) {
          if(data_values[j]==1){
            mini= max(mini,condi[j])
            
          }
          if(data_values[j]==0) {
            maxi=min(maxi,condi[j])
            
          }
        }
      }
      V = truncnorm(mini,maxi,0,1)
      sim = sim + c_data[i,]*(V-sim[i])
    }
  }
  for ( i in 1:n){
    x<-data_loc[i,1]
    y<-data_loc[i,2]
    z<-data_loc[i,3]
    grid[x,y,z]<-sim[i]
  }
  return(grid)
}

gibbs_sampling_diag<-function(grid,gridco,q,dz,dhx,dhy,rv,rhx,rhy,b,it){
  # this algorithm is based on Emery paper (2013); algorithm 7
  
  codata_values<-gridco[!is.na(gridco)] # getting the values
  data_values<-grid[!is.na(grid)]
  data_loc<-which(!is.na(grid),arr.ind=TRUE) # getting the position in the grid
  n=length(data_values)
  
  thres=matrix(0,2,n)
  thres[]=c(-Inf,Inf) # default, which works for the values 0 (well extended)
  thres[,codata_values==1 & data_values==2]=c(q[1],q[2])
  sima=truncnorm(q[1],q[2],0,1)
  thres[,codata_values==3 & data_values==2]=c(-Inf,q[1])
  simb=truncnorm(-Inf,q[1],0,1)
  thres[,codata_values==2 & data_values==2]=c(q[2],Inf)
  simc=truncnorm(q[2],Inf,0,1)
  thres[,codata_values==2 & data_values==3]=c(q[3],Inf)
  simd=truncnorm(q[3],Inf,0,1)
  thres[,codata_values==3 & data_values==3]=c(-Inf,q[3])
  sime=truncnorm(-Inf,q[3],0,1)
  # computation of the covariance matrix
  c_data<- cov_data(data_loc,dz,dhx,dhy,rhx,rhy,rv,b)
  #inv <- which(c_data <0)
  #initialization uncorrelated gaussian simulation inside the interval

  s=rnorm(1)
  Z<-rep(s,n)
  Z[codata_values==1 & data_values==2]=sima
  Z[codata_values==3 & data_values==2]=simb
  Z[codata_values==2 & data_values==2]=simc
  Z[codata_values==2 & data_values==3]=simd
  Z[codata_values==3 & data_values==3]=sime
  
  # gibbs sampling
  for( t in 1:it){
    for ( i in 1:n){
      condi_inf<-(((thres[1,]-Z)/c_data[i,])+Z[i])
      condi_sup<-(((thres[2,]-Z)/c_data[i,])+Z[i])
      mini=max(pmin(condi_inf,condi_sup),na.rm=TRUE) #eq(21) we can have NA when the covariance is 0 
      maxi=min(pmax(condi_inf,condi_sup),na.rm=TRUE)  #eq(22)
      #if (abs(pnorm(mini)-pnorm(maxi))<0.0001){print(c(mini,maxi,mean(Z),var(Z)))}
      #print(c(mini,maxi,mean(Z),var(Z)))
      V = truncnorm(mini,maxi,0,1)
      
      Z = Z + c_data[i,]*(V-Z[i])
    }
  }
  for ( i in 1:n){
    x<-data_loc[i,1]
    y<-data_loc[i,2]
    z<-data_loc[i,3]
    grid[x,y,z]<-Z[i]
  }
  return(grid)
}