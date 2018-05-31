# conditioning for four Gaussian functions
# the covariance matrix is going to be very big
# probably not inversible, so that we will probably need to use exponential covariances

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

coco_cov_data<-function(data_loc,data_loc_co,dz,dhx,rv,rhx,b,rho,shift){
  # the matrix is going to be non symmetric but positive definit
  # we are here coding the general case where n can be different from n_co
  n=nrow(data_loc)
  n_co=nrow(data_loc_co)
  dist_hx <- dist_v<- matrix(0,n,n)
  for ( i in 1:(n-1)){
    for (j in (i+1):n ){
      dist_v[i,j]=-(data_loc[i,2]-data_loc[j,2])*dz
      dist_v[j,i]=-dist_v[i,j]
      dist_hx[i,j]=dist_hx[j,i]=abs(data_loc[i,1]-data_loc[j,1])*dhx
    }
  }
  dist_hx_co <- dist_v_co<- matrix(0,n,n_co)
  co_dist_hx <- co_dist_hy <- co_dist_v<- matrix(0,n_co,n)
  for ( i in 1:n){
    for (j in 1:n_co ){
      dist_v_co[i,j]=-(data_loc[i,2]-data_loc_co[j,2])*dz
      co_dist_v[j,i]=-dist_v_co[i,j]
      dist_hx_co[i,j]=co_dist_hx[j,i]=abs(data_loc[i,1]-data_loc_co[j,1])*dhx
    }
  }
  co_dist_hx_co  <- co_dist_v_co<- matrix(0,n_co,n_co)
  for ( i in 1:(n_co-1)){
    for (j in (i+1):n_co ){
      co_dist_v_co[i,j]=-(data_loc_co[i,2]-data_loc_co[j,2])*dz
      co_dist_v_co[j,i]=-co_dist_v_co[i,j]
      co_dist_hx_co[i,j]=co_dist_hx_co[j,i]=abs(data_loc_co[i,1]-data_loc_co[j,1])*dhx
    }
  }
  a12= rho[1]/model(shift[1],rv[1],b[1])
  a14=rho[2]/model(shift[2],rv[1],b[1])
  a24=rho[3]/model(shift[3],rv[2],b[2])
  #horizontal covariances
  cov11x=model_lat(dist_hx,rhx[1])
  cov11x_co=model_lat(dist_hx_co,rhx[1])
  co_cov11x=model_lat(co_dist_hx,rhx[1])

  
  cov12=a12*model(dist_v+shift[1],rv[1],b[1])*cov11x
  cov14=a14*model(dist_v_co+shift[2],rv[1],b[1])*cov11x_co
  
  cov21=a12*model(dist_v-shift[1],rv[1],b[1])*cov11x
  cov41=a14*model(co_dist_v-shift[2],rv[1],b[1])*co_cov11x
  
  cov24 = ( a12*a14*model(dist_v_co-shift[1]+shift[2],rv[1],b[1]) + sqrt(1-a12^2)*a24*model(dist_v_co+shift[3],rv[2],b[2]) )*cov11x_co
  cov42 = (a12*a14*model(co_dist_v+shift[1]-shift[2],rv[1],b[1])+ sqrt(1-a12^2)*a24*model(co_dist_v-shift[3],rv[2],b[2]) )*co_cov11x
  
  cov11=model(dist_v,rv[1],b[1])*cov11x
  cov22=(a12^2)*model(dist_v,rv[1],b[1])*cov11x+(1-(a12)^2)*model_lat(dist_hx,rhx[2])*model(dist_v,rv[2],b[2])
  cov33=model(co_dist_v_co,rv[3],b[3])*model_lat(co_dist_hx_co,rhx[3])
  cov44=(a14^2)*model(co_dist_v_co,rv[1],b[1])*model_lat(co_dist_hx_co,rhx[1])+(a24^2)*model(co_dist_v_co,rv[2],b[2])*model_lat(co_dist_hx_co,rhx[2])+(1-(a14^2)-(a24^2))*model_lat(co_dist_hx_co,rhx[4])*model(co_dist_v_co,rv[4],b[4])
  
  # now the cross correlations that equals 0
  cov34=cov13=cov23=0*dist_v_co
  cov43=cov32=cov31=0*co_dist_v
  
  #building the final matrix as a concatenation! 
  return(rbind(cbind(cov11,cov12,cov13,cov14),cbind(cov21,cov22,cov23,cov24),cbind(cov31,cov32,cov33,cov34),cbind(cov41,cov42,cov43,cov44)))
}

coco_gibbs_sampling<-function(grid,grid_co,q,q_co,dz,dhx,rv,rhx,b,rho,shift,it){
  # this algorithm is based on Emery paper (2013); algorithm 7
  
  data_values<-grid[!is.na(grid)] # getting the values
  data_values_co<-grid_co[!is.na(grid_co)]
  data_loc<-which(!is.na(grid),arr.ind=TRUE) # getting the position in the grid
  data_loc_co<-which(!is.na(grid_co),arr.ind=TRUE)
  n=length(data_values)
  n_co=length(data_values_co) # we have to be careful, because the number of points for the diagenesis might not be the same!
  
  # defining the thresholds for each Gaussian function
  thres1=thres2=matrix(0,2,n)
  thres1[]=thres2[]=c(-Inf,Inf) # default, which works for the values 0 (well extended)

  thres1[,data_values==1]=c(-Inf,q[1])
  thres2[,data_values==1]=c(-Inf,q[2])
  
  thres1[,data_values==2]=c(-Inf,q[1])
  thres2[,data_values==2]=c(q[2],Inf)
  
  thres1[,data_values==3]=c(q[1],Inf)
  thres2[,data_values==3]=c(q[3],Inf)
  
  thres1[,data_values==4]=c(q[1],Inf)
  thres2[,data_values==4]=c(-Inf,q[3])
  
  thres=cbind(thres1,thres2)
  thres3=thres4=matrix(0,2,n_co)
  thres3[]=thres4[]=c(-Inf,Inf) # default
  thres4[,data_values_co==3]=c(q_co[1],Inf)
  thres4[,data_values_co==2 & data_values==2]=c(q_co[2],Inf)
  thres4[,data_values_co==2 & data_values==3]=c(q_co[3],Inf)
  
  thres4[,data_values_co==1 & data_values==2]=c(-Inf,q_co[2])
  thres4[,data_values_co==1 & data_values==3]=c(-Inf,q_co[3])
  thres4[,data_values_co==1 & data_values==4]=c(-Inf,q_co[1])

  thres=cbind(thres,thres3,thres4)
  # computation of the co-kriging covariance matrix, with four variables !! 
  c_data<- coco_cov_data(data_loc,data_loc_co,dz,dhx,rv,rhx,b,rho,shift)
  
  #initialization the two gaussian functions (uncorrelated to each other)
  Z<-rep(0,2*n+2*n_co)
  for ( i in 1:n){
    if (data_values[i]==1){
      Z[i] =  negrnorm(q[1],0,1)
      Z[i+n] = negrnorm(q[2],0,1)
    }
    if (data_values[i]==2){
      Z[i] = negrnorm(q[1],0,1)
      Z[i+n]  = posrnorm(q[2],0,1)
    }
    if (data_values[i]==3){
      Z[i] = posrnorm(q[1],0,1)
      Z[i+n]  = posrnorm(q[3],0,1)
    }
    if (data_values[i]==4){
      Z[i] = posrnorm(q[1],0,1)
      Z[i+n] = negrnorm(q[3],0,1)
    }
    if (data_values[i]==0){ # no boundaries
      Z[i] =  rnorm(1)
      Z[i+n] = rnorm(1)
    }
    if (data_values_co[i]==3){
      Z[i+2*n] = rnorm(1) 
      Z[i+2*n+n_co] = posrnorm(q_co[1],0,1)
    }
    if (data_values_co[i]==2 & data_values[i]==2 ){
      Z[i+2*n] =  rnorm(1)
      Z[i+2*n+n_co] = posrnorm(q_co[2],0,1)
    }
    if (data_values_co[i]==2 & data_values[i]==3 ){
      Z[i+2*n] =  rnorm(1)
      Z[i+2*n+n_co] = posrnorm(q_co[3],0,1)
    }
    if (data_values_co[i]==1 & data_values[i]==2 ){
      Z[i+2*n] =  rnorm(1)
      Z[i+2*n+n_co] = negrnorm(q_co[2],0,1)
    }
    if (data_values_co[i]==1 & data_values[i]==3 ){
      Z[i+2*n] =  rnorm(1)
      Z[i+2*n+n_co] = negrnorm(q_co[3],0,1)
    }
    if (data_values_co[i]==1 & data_values[i]==4 ){
      Z[i+2*n] =  negrnorm(q_co[1],0,1)
      Z[i+2*n+n_co] = rnorm(1)
    }
    if (data_values_co[i]==1 & data_values[i]==1 ){
      Z[i+2*n] =  rnorm(1)
      Z[i+2*n+n_co] = rnorm(1)
    }
    if (data_values_co[i]==0){ # no boundaries
      Z[i+2*n] =  rnorm(1)
      Z[i+2*n+n_co] = rnorm(1)
    }
  }
  
  # gibbs sampling
  for( t in 1:it){
    for ( i in 1:(2*n+2*n_co)){
      # equations from emery et al. (2014)
      condi_inf<-(((thres[1,]-Z)/c_data[i,])+Z[i])
      condi_sup<-(((thres[2,]-Z)/c_data[i,])+Z[i])
      mini=max(pmin(condi_inf,condi_sup),na.rm=TRUE) #eq(21) we can have NA when the covariance is 0 
      maxi=min(pmax(condi_inf,condi_sup),na.rm=TRUE)  #eq(22)
      #if (abs(pnorm(mini)-pnorm(maxi))<0.0001){print(c(mini,maxi,mean(Z),var(Z)))}
      #print(c(mini,maxi,mean(Z),var(Z)))
      V = truncnorm(mini,maxi,0,1)
      
      Z = Z + c_data[i,]*(V-Z[i]) # this algorithm allows to modify all the values of the simulation at every step 
    }
  }
  Z1<-Z2<-Z3<-Z4<-grid
  for ( i in 1:n){
    z<-data_loc[i,2]
    x<-data_loc[i,1]
    z_co<-data_loc_co[i,2]
    x_co<-data_loc_co[i,1]
    Z1[x,z]<-Z[i]
    Z2[x,z]<-Z[i+n]
    Z3[x_co,z_co]<-Z[i+2*n]
    Z4[x_co,z_co]<-Z[i+2*n+n_co]
  }
  return(cbind(Z1,Z2,Z3,Z4))
}

coco_dual_krig<-function(E1well,E2well,E3well,E4well,dz,dhx,rv,rhx,b,rho,shift){
  
  #positions of data and estimate points
  data_loc<-which(!is.na(E1well),arr.ind=TRUE) #we make the assumption here that all wells have same position, which is not always true
  estim_loc <-which(is.na(E1well),arr.ind=TRUE)
  data_values=c(E1well[!is.na(E1well)],E2well[!is.na(E2well)],E3well[!is.na(E3well)],E4well[!is.na(E4well)])
  m = nrow(estim_loc)
  
  c_data<- coco_cov_data(data_loc,data_loc,dz,dhx,rv,rhx,b,rho,shift)
  s<-solve(c_data,data_values) # it's better to solve the system rather than having the matrix in memory (and more accurate)
  #is.positive.definite(c_data)
  a12= rho[1]/model(shift[1],rv[1],b[1])
  a14=rho[2]/model(shift[2],rv[1],b[1])
  a24=rho[3]/model(shift[3],rv[2],b[2])
  
  for (i in 1:m){
    xi=estim_loc[i,1]
    zi=estim_loc[i,2]
    
    dist_hx=abs(data_loc[,1]-xi)*dhx
    dist_v= -(data_loc[,2]-zi)*dz # the point to estimate is like x+h, the sign here is quite important
    
    #horizontal covariances
    cov11x=model_lat(dist_hx,rhx[1])
    
    cov12=a12*model(dist_v+shift[1],rv[1],b[1])*cov11x
    cov14=a14*model(dist_v+shift[2],rv[1],b[1])*cov11x
    
    cov21=a12*model(dist_v-shift[1],rv[1],b[1])*cov11x
    cov41=a14*model(dist_v-shift[2],rv[1],b[1])*cov11x
    
    cov24 = ( a12*a14*model(dist_v-shift[1]+shift[2],rv[1],b[1]) + sqrt(1-a12^2)*a24*model(dist_v+shift[3],rv[2],b[2]) )*cov11x
    cov42 = (a12*a14*model(dist_v+shift[1]-shift[2],rv[1],b[1])+ sqrt(1-a12^2)*a24*model(dist_v-shift[3],rv[2],b[2]) )*cov11x
    
    cov11=model(dist_v,rv[1],b[1])*cov11x
    cov22=(a12^2)*model(dist_v,rv[1],b[1])*cov11x+(1-(a12)^2)*model_lat(dist_hx,rhx[2])*model(dist_v,rv[2],b[2])
    cov33=model(dist_v,rv[3],b[3])*model_lat(dist_hx,rhx[3])
    cov44=(a14^2)*model(dist_v,rv[1],b[1])*cov11x
    +(a24^2)*model(dist_v,rv[2],b[2])*model_lat(dist_hx,rhx[2])
    +(1-(a14^2)-(a24^2))*model_lat(dist_hx,rhx[4])*model(dist_v,rv[4],b[4])
    
    # now the cross correlations that equals 0
    cov34=cov13=cov23=cov43=cov32=cov31=0.*dist_v
    
    E1well[xi,zi]= s %*% c(cov11,cov12,cov13,cov14)
    E2well[xi,zi]= s %*% c(cov21,cov22,cov23,cov24)
    E3well[xi,zi]= s %*% c(cov31,cov32,cov33,cov34)
    E4well[xi,zi]= s %*% c(cov41,cov42,cov43,cov44)
  }
  
  
  return(cbind(E1well,E2well,E3well,E4well))
}


cov_data_surf<-function(data_loc,dhx,dhy,rhx,rhy){
  n=nrow(data_loc)
  dist_hx <-dist_hy <- matrix(0,n,n)
  for ( i in 1:(n-1)){
    for (j in (i+1):n ){
      dist_hx[i,j]=dist_hx[j,i]=abs(data_loc[i,1]-data_loc[j,1])*dhx
      dist_hy[i,j]=dist_hy[j,i]=abs(data_loc[i,2]-data_loc[j,2])*dhy
    }
  }
  
  cov_h=exp(-sqrt( ((dist_hx/rhx)^2) + ((dist_hy/rhy)^2) ))
  
  
  return(cov_h)
}

cov_data_line<-function(data_loc,dhx,rhx){
  n=length(data_loc)
  dist_hx <-dist_hy <- matrix(0,n,n)
  for ( i in 1:(n-1)){
    for (j in (i+1):n ){
      dist_hx[i,j]=dist_hx[j,i]=abs(data_loc[i]-data_loc[j])*dhx
    }
  }
  
  cov_h=exp(-sqrt( ((dist_hx/rhx)^2) ))
  
  
  return(cov_h)
}


dual_krig_surf<-function(errorwell,rhx,rhy,dhx,dhy){
  
  # getting the positions in the grid of one horizontal line, the positions will be the same for every line
  data_loc<-which(!is.na(errorwell[,,1]),arr.ind=TRUE) 
  toestim <-which(is.na(errorwell[,,1]),arr.ind=TRUE)
  n= nrow(data_loc)
  m = nrow(toestim)
  print(n)
  # covariance and inverse covariance matrix of a single surface
  c_data<- cov_data_surf(data_loc,dhx,dhy,rhx,rhy)
  invc<-chol2inv(chol(c_data))
  errorfield <-errorwell
  print(c_data)
  print(invc)
  
  # covariance vector between estimation and data of every surface
  c_estim_data<-matrix(0,m,n)
  for ( j in 1:m){
    hx=(data_loc[,1]-toestim[j,1])*dhx
    hy=(data_loc[,2]-toestim[j,2])*dhy
    c_estim_data[j,]<-exp(-sqrt( ((hx/rhx)^2)+ ((hy/rhy)^2) ))
  }
  
  # for every surface
  ndim<-dim(errorwell)
  for ( i in 1:ndim[3]){
    surf<-errorwell[,,i]
    data_values<-surf[!is.na(surf)]
    b <- t(data_values) %*% invc
    for ( j in 1:m){
      hx = toestim[j,1]
      hy = toestim[j,2]
      errorfield[hx,hy,i] = b %*% c_estim_data[j,]
    }
  }
  
  
  return(errorfield)
}

dual_krig_line<-function(errorwell,rhx,dhx){
  
  # getting the positions in the grid of one horizontal line, the positions will be the same for every line
  data_loc<-which(!is.na(errorwell[,1]),arr.ind=TRUE) 
  toestim <-which(is.na(errorwell[,1]),arr.ind=TRUE)
  n= length(data_loc)
  m = length(toestim)
  print(n)
  # covariance and inverse covariance matrix of a single surface
  c_data<- cov_data_line(data_loc,dhx,rhx)
  invc<-chol2inv(chol(c_data))
  errorfield <-errorwell
  print(c_data)
  print(invc)
  
  # covariance vector between estimation and data of every surface
  c_estim_data<-matrix(0,m,n)
  for ( j in 1:m){
    hx=(data_loc[]-toestim[j])*dhx
    c_estim_data[j,]<-exp(-sqrt( ((hx/rhx)^2) ))
  }
  
  # for every surface
  for ( i in 1:ncol(errorwell)){
    surf<-errorwell[,i]
    data_values<-surf[!is.na(surf)]
    b <- t(data_values) %*% invc
    for ( j in 1:m){
      hx = toestim[j]
      errorfield[hx,i] = b %*% c_estim_data[j,]
    }
  }
  
  
  return(errorfield)
}


gibbs_sampling3d<-function(grid,grid_co,q,q_co,dz,dhx,dhy,rv,rhx,rhy,b,rho,shift,it){
  # this algorithm is based on Emery paper (2013); algorithm 7
  
  data_values<-grid[!is.na(grid)] # getting the values
  data_values_co<-grid_co[!is.na(grid_co)]
  data_loc<-which(!is.na(grid),arr.ind=TRUE) # getting the position in the grid
  data_loc_co<-which(!is.na(grid_co),arr.ind=TRUE)
  n=length(data_values)
  n_co=length(data_values_co) # we have to be careful, because the number of points for the diagenesis might not be the same!
  
  # defining the thresholds for each Gaussian function
  thres1=thres2=matrix(0,2,n)
  thres1[]=thres2[]=c(-Inf,Inf) # default, which works for the values 0 (well extended)
  
  thres1[,data_values==1]=c(-Inf,q[1])
  thres2[,data_values==1]=c(-Inf,q[2])
  
  thres1[,data_values==2]=c(-Inf,q[1])
  thres2[,data_values==2]=c(q[2],Inf)
  
  thres1[,data_values==3]=c(q[1],Inf)
  thres2[,data_values==3]=c(q[3],Inf)
  
  thres1[,data_values==4]=c(q[1],Inf)
  thres2[,data_values==4]=c(-Inf,q[3])
  
  thres=cbind(thres1,thres2)
  thres3=thres4=matrix(0,2,n_co)
  thres3[]=thres4[]=c(-Inf,Inf) # default
  thres4[,data_values_co==3]=c(q_co[1],Inf)
  thres4[,data_values_co==2 & data_values==2]=c(q_co[2],Inf)
  thres4[,data_values_co==2 & data_values==3]=c(q_co[3],Inf)
  
  thres4[,data_values_co==1 & data_values==2]=c(-Inf,q_co[2])
  thres4[,data_values_co==1 & data_values==3]=c(-Inf,q_co[3])
  thres4[,data_values_co==1 & data_values==4]=c(-Inf,q_co[1])
  
  thres=cbind(thres,thres3,thres4)
  # computation of the co-kriging covariance matrix, with four variables !! 
  c_data<- coco_cov_data3d(data_loc,data_loc_co,dz,dhx,dhy,rv,rhx,rhy,b,rho,shift)
  
  #initialization the two gaussian functions (uncorrelated to each other)
  # it is better to use one value for each interval
  za=negrnorm(q[1],0,1)
  zb=negrnorm(q[2],0,1)
  zc=posrnorm(q[2],0,1)
  zd=posrnorm(q[1],0,1)
  ze=posrnorm(q[3],0,1)
  zf=negrnorm(q[3],0,1)
  zg=rnorm(1)
  Z<-rep(0,2*n+2*n_co)
  for ( i in 1:n){
    if (data_values[i]==1){
      Z[i] =  za
      Z[i+n] = zb
    }
    if (data_values[i]==2){
      Z[i] = za
      Z[i+n]  = zc
    }
    if (data_values[i]==3){
      Z[i] = zd
      Z[i+n]  = ze
    }
    if (data_values[i]==4){
      Z[i] = zd
      Z[i+n] = zf
    }
    if (data_values[i]==0){ # no boundaries
      Z[i] =  zg
      Z[i+n] = zg
    }
    # we need to be careful now because only three gaussian functions are used
    # and the fourth is used 
    zh=posrnorm(q_co[1],0,1)
    zi=posrnorm(q_co[2],0,1)
    zj=posrnorm(q_co[3],0,1)
    zk=negrnorm(q_co[2],0,1)
    zl=negrnorm(q_co[3],0,1)
    zm=negrnorm(q_co[1],0,1)
    if (data_values_co[i]==3){
      Z[i+2*n] =  zg
      Z[i+2*n+n_co] = zh
    }
    if (data_values_co[i]==2 & data_values[i]==2 ){
      Z[i+2*n] =  zg
      Z[i+2*n+n_co] = zi
    }
    if (data_values_co[i]==2 & data_values[i]==3 ){
      Z[i+2*n] =  zg
      Z[i+2*n+n_co] = zj
    }
    if (data_values_co[i]==1 & data_values[i]==2 ){
      Z[i+2*n] =  zg
      Z[i+2*n+n_co] = zk
    }
    if (data_values_co[i]==1 & data_values[i]==3 ){
      Z[i+2*n] =  zg
      Z[i+2*n+n_co] = zl
    }
    if (data_values_co[i]==1 & data_values[i]==4 ){
      Z[i+2*n] =  zg
      Z[i+2*n+n_co] = zm
    }
    if (data_values_co[i]==1 & data_values[i]==1 ){
      Z[i+2*n] =  zg
      Z[i+2*n+n_co] = zg
    }
    if (data_values_co[i]==0){ # no boundaries
      Z[i+2*n] =  zg
      Z[i+2*n+n_co] = zg
    }
  }
  
  # gibbs sampling
  for( t in 1:it){
    for ( i in 1:(2*n+2*n_co)){
      # equations from emery et al. (2014)
      condi_inf<-(((thres[1,]-Z)/c_data[i,])+Z[i])
      condi_sup<-(((thres[2,]-Z)/c_data[i,])+Z[i])
      mini=max(pmin(condi_inf,condi_sup),na.rm=TRUE) #eq(21) we can have NA when the covariance is 0 
      maxi=min(pmax(condi_inf,condi_sup),na.rm=TRUE)  #eq(22)
      #if (abs(pnorm(mini)-pnorm(maxi))<0.0001){print(c(mini,maxi,mean(Z),var(Z)))}
      #print(c(mini,maxi,mean(Z),var(Z)))
      V = truncnorm(mini,maxi,0,1)
      #print(c(mini,maxi))
      
      Z = Z + c_data[i,]*(V-Z[i]) # this algorithm allows to modify all the values of the simulation at every step 
    }
  }
  Z1<-Z2<-Z3<-Z4<-grid
  for ( i in 1:n){
    z<-data_loc[i,3]
    y<-data_loc[i,2]
    x<-data_loc[i,1]
    z_co<-data_loc_co[i,3]
    y_co<-data_loc_co[i,2]
    x_co<-data_loc_co[i,1]
    Z1[x,y,z]<-Z[i]
    Z2[x,y,z]<-Z[i+n]
    Z3[x_co,y_co,z_co]<-Z[i+2*n]
    Z4[x_co,y_co,z_co]<-Z[i+2*n+n_co]
  }
  return(abind(Z1,Z2,Z3,Z4,along=3))
}

coco_cov_data3d<-function(data_loc,data_loc_co,dz,dhx,dhy,rv,rhx,rhy,b,rho,shift){
  # the matrix is going to be non symmetric but positive definit
  # we are here coding the general case where n can be different from n_co
  n=nrow(data_loc)
  n_co=nrow(data_loc_co)
  dist_hx <-dist_hy<- dist_v<- matrix(0,n,n)
  for ( i in 1:(n-1)){
    for (j in (i+1):n ){
      dist_v[i,j]=-(data_loc[i,3]-data_loc[j,3])*dz
      dist_v[j,i]=-dist_v[i,j]
      dist_hx[i,j]=dist_hx[j,i]=abs(data_loc[i,1]-data_loc[j,1])*dhx
      dist_hy[i,j]=dist_hy[j,i]=abs(data_loc[i,2]-data_loc[j,2])*dhy
    }
  }
  dist_hx_co<-dist_hy_co <- dist_v_co<- matrix(0,n,n_co)
  co_dist_hx <- co_dist_hy <- co_dist_v<- matrix(0,n_co,n)
  for ( i in 1:n){
    for (j in 1:n_co ){
      dist_v_co[i,j]=-(data_loc[i,3]-data_loc_co[j,3])*dz
      co_dist_v[j,i]=-dist_v_co[i,j]
      dist_hx_co[i,j]=co_dist_hx[j,i]=abs(data_loc[i,1]-data_loc_co[j,1])*dhx
      dist_hy_co[i,j]=co_dist_hy[j,i]=abs(data_loc[i,2]-data_loc_co[j,2])*dhy
    }
  }
  co_dist_hx_co <-co_dist_hy_co <- co_dist_v_co<- matrix(0,n_co,n_co)
  for ( i in 1:(n_co-1)){
    for (j in (i+1):n_co ){
      co_dist_v_co[i,j]=-(data_loc_co[i,3]-data_loc_co[j,3])*dz
      co_dist_v_co[j,i]=-co_dist_v_co[i,j]
      co_dist_hx_co[i,j]=co_dist_hx_co[j,i]=abs(data_loc_co[i,1]-data_loc_co[j,1])*dhx
      co_dist_hy_co[i,j]=co_dist_hy_co[j,i]=abs(data_loc_co[i,2]-data_loc_co[j,2])*dhy
    }
  }
  a12= rho[1]/model(shift[1],rv[1],b[1])
  a14=rho[2]/model(shift[2],rv[1],b[1])
  a24=rho[3]/model(shift[3],rv[2],b[2])
  #horizontal covariances
  cov11x=model_lat(dist_hx,rhx[1])
  cov11y=model_lat(dist_hy,rhy[1])
  cov11x_co=model_lat(dist_hx_co,rhx[1])
  cov11y_co=model_lat(dist_hy_co,rhy[1])
  co_cov11x=model_lat(co_dist_hx,rhx[1])
  co_cov11y=model_lat(co_dist_hy,rhy[1])
  
  
  cov12=a12*model(dist_v+shift[1],rv[1],b[1])*cov11x*cov11y # this one is separable
  cov14=a14*model(dist_v_co+shift[2],rv[1],b[1])*cov11x_co*cov11y_co
  
  cov21=a12*model(dist_v-shift[1],rv[1],b[1])*cov11x*cov11y
  cov41=a14*model(co_dist_v-shift[2],rv[1],b[1])*co_cov11x*co_cov11y
  
  cov24 = ( a12*a14*model(dist_v_co-shift[1]+shift[2],rv[1],b[1]) + sqrt(1-a12^2)*a24*model(dist_v_co+shift[3],rv[2],b[2]) )*cov11x_co*cov11y_co
  cov42 = (a12*a14*model(co_dist_v+shift[1]-shift[2],rv[1],b[1])+ sqrt(1-a12^2)*a24*model(co_dist_v-shift[3],rv[2],b[2]) )*co_cov11x*co_cov11y
  
  cov11=model(dist_v,rv[1],b[1])*cov11x*cov11y
  cov22=(a12^2)*model(dist_v,rv[1],b[1])*cov11x*cov11y+(1-(a12)^2)*model_lat(dist_hx,rhx[2])*model_lat(dist_hy,rhy[2])*model(dist_v,rv[2],b[2])
  cov33=model(co_dist_v_co,rv[3],b[3])*model_lat(co_dist_hx_co,rhx[3])*model_lat(co_dist_hy_co,rhy[3])
  cov44=(a14^2)*model(co_dist_v_co,rv[1],b[1])*model_lat(co_dist_hx_co,rhx[1])*model_lat(co_dist_hy_co,rhy[1])+(a24^2)*model(co_dist_v_co,rv[2],b[2])*model_lat(co_dist_hx_co,rhx[2])*model_lat(co_dist_hy_co,rhy[2])+(1-(a14^2)-(a24^2))*model_lat(co_dist_hx_co,rhx[4])*model_lat(co_dist_hy_co,rhy[4])*model(co_dist_v_co,rv[4],b[4])
  
  # now the cross correlations that equals 0
  cov34=cov13=cov23=0*dist_v_co
  cov43=cov32=cov31=0*co_dist_v
  
  #building the final matrix as a concatenation! 
  return(rbind(cbind(cov11,cov12,cov13,cov14),cbind(cov21,cov22,cov23,cov24),cbind(cov31,cov32,cov33,cov34),cbind(cov41,cov42,cov43,cov44)))
}