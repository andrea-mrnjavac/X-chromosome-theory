NeA<-100
NeX<-75

MV_A<-function(x,h,sm,sf){
  return(NeA*(sm+sf)*(x*(1-2*h)+h)) 
}

MV_hX<-function(x,h,sm,sf){
  return((2/3)*NeX*(sf*(2*h+2*x*(1-2*h))+sm))}

MV_dX<-function(x,h,sm,sf){
  return((2/3)*NeX*(sf*(2*h+2*x*(1-2*h))+ h*sm)) }

#autosomal loci

# function G
G_A<-function(xi,h,sm,sf){
  int<-integrate(MV_A,0,xi,h=h,sm=sm, sf=sf)
  gx<-exp(-2*int$value)
  return(gx)
}

#function U
U_A<-function(p,h,sm,sf){
  int1<-integrate(Vectorize(G_A),0,p,h=h,sm=sm,sf=sf)
  int2<-integrate(Vectorize(G_A),0,1,h=h,sm=sm,sf=sf)
  u<-int1$value/int2$value
  return(u)
}

#hemizygous X loci

G_h<-function(xi,h,sm,sf){
  int<-integrate(MV_hX,0,xi,h=h,sm=sm, sf=sf)
  gx<-exp(-2*int$value)
  return(gx)
}

U_h<-function(p,h,sm,sf){
  int1<-integrate(Vectorize(G_h),0,p,h=h,sm=sm,sf=sf)
  int2<-integrate(Vectorize(G_h),0,1,h=h,sm=sm,sf=sf)
  u<-int1$value/int2$value
  return(u)
}

#diploid X loci

G_d<-function(xi,h,sm,sf){
  int<-integrate(MV_dX,0,xi,h=h,sm=sm, sf=sf)
  gx<-exp(-2*int$value)
  return(gx)
}

U_d<-function(p,h,sm,sf){
  int1<-integrate(Vectorize(G_d),0,p,h=h,sm=sm,sf=sf)
  int2<-integrate(Vectorize(G_d),0,1,h=h,sm=sm,sf=sf)
  u<-int1$value/int2$value
  return(u)
}
