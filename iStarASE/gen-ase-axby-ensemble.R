# Compute an ensemble with members that have fields like
# topg <- topg_mb(x,y) + f_topg(x,y)*sigma_topg(x,y)
# thck <- thck_mb(x,y) + f_thck(x,y)*sigma_thck(x,y)
# C <- C_mb * exp(f_C(x,y) * sigma_lnC(x,y)) [basal friction coefficient]
# D <- D_mb * exp(f_D(x,y) * sigma_lnD(x,y)) [viscosity factor]
# M <- M_mb * exp(f_D(x,y) * sigma_lnD(x,y)) [melt rate]
amr.free.all()
#eneumerate components 

TOPG <- 0
THCK <-  1
CWSHELF <- 2
C <- 3
MUCOEF <- 4
MELT <- 5

#the f_i(x,y) will be contstructed by adding hat functions on chosen centers
#(e.g the two parts of Thwaites shelf)

#the sigma_i are constructed either from the difference between the
#inverse problem (with mass balance) topography/thickness (topg_mb/thick_mb)
#and the bedmap2 topgraphy/thickness data, or by choosing sane values


require(libamrfile)
require(lhs)

amr.aXbYplusc <- function(X,Y,a,b=function(x,y,comp){1-a(x,y,comp)},c=function(x,y,comp){0})
{

  #given two AMR files X and Y , produce an amr hierarchy one whose components are
  #a(x,y,comp) * X + a(b,y,comp) * Y + c(x,y,comp) 
  
  #ASSUMES X,Y HAVE THE SAME DISJOINT BOX LAYOUT!
   xID <- amr.load(X)
   yID <- amr.load(Y)

   maxlev <- amr.query.nlevel(xID) - 1
   maxlevy <- amr.query.nlevel(yID) - 1
   if (maxlev != maxlevy)
     {
       stop ("amr.axby assumes identical meshes but level count differs")
     }

   for (lev in 0:maxlev)
     {
       nfab <- amr.query.nfab(xID,lev) - 1
       nfaby <- amr.query.nfab(yID,lev) - 1
       if (nfab != nfaby)
         {
           stop ("amr.axby assumes identical meshes but fab count differs")
         }

       ncomp <- amr.query.ncomp(xID, lev)
       ncompy <- amr.query.ncomp(yID, lev)
       if (ncomp != ncompy)
         {
           stop ("amr.axby assumes same number of components")
         }

       for (ifab in 0:nfab)
         {
           for (icomp in 0:(ncomp-1))
             {
               fab <- amr.read.fab(xID,lev,ifab,icomp,ng=1)
               faby <- amr.read.fab(yID,lev,ifab,icomp,ng=1)

               data <- a(fab$x,fab$y,icomp) * fab$v +
                 b(fab$x,fab$y,icomp) * faby$v +
                   c(fab$x,fab$y,icomp)
               
               
               amr.write.fab(xID, lev, ifab, data,icomp , ng=1)
               
             }
         }
     }
   
   print(xID)
   ncomp <- amr.query.ncomp(xID, 0)
   print(ncomp)
   print(yID)
   ncomp <- amr.query.ncomp(xID, 0)
   print(ncomp)
   amr.free(yID)
   print(xID)
   ncomp <- amr.query.ncomp(xID, 0)
   print(ncomp)
   xID
 }


create_sigma <- function(){
#create the sigma data for thickness and topography
#(MAYBE WE SHOULD STORE TOPG,THCK sperately...)
  
#contruct a,b,c such that we can use amr.aXbYplusc to find
#sigma_topg =  topg_bedmap - topg_bm
#sigma_thck =  thck_bedmap - thck_bm
#sigma_lnC = 0.0 (MAYBE WE SHOULD STORE TOPG,THCK sperately...)
#sigma_lnD = 0.0
#sigma_M =   0.0

  a <- function(x,y,comp)
    {
      if (comp ==THCK | comp == TOPG ){
                                        #thickness and topography data
        1.0
      } else {
        0.0
      } 
  }
  b <- function(x,y,comp){-a(x,y,comp)}
  
  c <- function(x,y,comp)
    {
      matrix(0.0,length(x),length(y))
   }
  
    
  amrID <- amr.aXbYplusc("ase-data-mb-1km.2d.hdf5","ase-data-zero-1km.2d.hdf5", a,b,c )
  ncomp <- amr.query.ncomp(amrID, lev)
  amr.write(amrID, "ase-data-sigma-1km.2d.hdf5")
  amr.free(amrID)
}



hat <- function (x,y,origin = c(0.0,0.0), width = 1.0, a = 1.0)
{
   #2D hat function h(x,y)
   #h(x,y) = a at the origin, zero for x > origin[0] +- width, etc
  print(origin); print(width); print (range(y))
   x <- (x - origin[1])/(width)
   y <- (y - origin[2])/(width)
   nx <- length(x)
   ny <- length(y)
   
   xx <- abs(matrix(rep(x,ny),nx,ny))
   yy <- abs(matrix(rep(y,each=nx),nx,ny))
   
   xk <- ifelse(xx < 1, 1 - abs(xx), 0)
   yk <- ifelse(yy < 1, 1 - abs(yy), 0)

   a * xk * yk
}

haty <- function (x,y,origin = c(0.0,0.0), width = 1.0, a = 1.0)
{
   #2D hat-in-y function h(x,y)
   #h(x,y) = a at the origin, zero for y > origin[0] +- width, etc
  print(origin); print(width); print (range(y))
   x <- (x - origin[1])/(width)
   y <- (y - origin[2])/(width)
   nx <- length(x)
   ny <- length(y)
   
   xx <- abs(matrix(rep(x,ny),nx,ny))
   yy <- abs(matrix(rep(y,each=nx),nx,ny))
   
   xk <- ifelse(xx < 1, 1 , 1)
   yk <- ifelse(yy < 1, 1 - abs(yy), 0)

   a * xk * yk
}

thwaites_pert_factory <- function(p)
  {
    #map vector p to p/2 functions, each of which is
    #a weighted sum over hat functions centered on different
    #parts of the thwaites gl
    # distance between nodes 
    hatwidth <- 64e+3 
    #centers of the perturbed regions
    thwaitesE <- c(4.5,7.5)*hatwidth
    thwaitesW = c(4.5,6.5)*hatwidth

    #sum of two hats
    g <- function(x,y,a,b)
      {
        
        a * hat(x,y,origin=thwaitesE,width=hatwidth) +
          b * hat(x,y,origin=thwaitesW,width=hatwidth) 
      }

      #sum of two y-hats
    gy <- function(x,y,a,b)
      {
        
        a * haty(x,y,origin=thwaitesE,width=hatwidth) +
          b * haty(x,y,origin=thwaitesW,width=hatwidth) 
      }
    
    list( Hf = function(x,y){gy(x,y,p[1],p[2])},
          Cf = function(x,y){gy(x,y,p[3],p[4])},
          Df = function(x,y){gy(x,y,p[5],p[6])},
          Mf = function(x,y){gy(x,y,p[7],p[8])})
  }


apply_thawites_pert <- function(p, X, Y)
  {

    if (length(p) != 8)
      {
        stop ("number of parameters must be 8 : this function modifies 4 fields in 2 location")
      }

    f <- thwaites_pert_factory(p)

    #construct aXbYplusc functions
    af <- function(x,y,comp)
      {
        r <- 1
        if (comp == CWSHELF | comp == C) 
          {
            #basal traction multiplier
            r <- 2^(f$Cf(x,y))
          }
        if (comp == MUCOEF) 
          {
            #viscosity multiplier
            r <- 2^(f$Df(x,y))
          }
        r
      }
    
    bf <- function(x,y,comp)
      {
         r <- 0
         
        if (comp == THCK | comp == TOPG)
          {
            #thickness or topography : add hat * Y
            r <-  f$Hf(x,y)
          }
        r
      }
              
    cf <- function(x,y,comp)
      {
        r <- matrix(0.0,length(x),length(y))
        if (comp == MELT)
          {
            r <- 100    * f$Mf(x,y)
            
          }
        r
      }

    amr.aXbYplusc("ase-data-mb-1km.2d.hdf5","ase-data-sigma-1km.2d.hdf5", af,bf,cf )
    
  }


if (FALSE)
  {
    create_sigma()
  }

p <- rep(c(1,-1),4)
f <- thwaites_pert_factory(p)

x <- seq(0,896e3,by=1e4)
y <- seq(0,1024e3,by=1e4)

par(mfrow=c(2,2))

#image(f$Hf(x,y))
#image(f$Df(x,y))
#image(f$Mf(x,y))
#image(f$Cf(x,y))

pertID <- apply_thawites_pert(p)
amr.write(pertID, "ase-data-pert.2d.hdf5")
amr.free(pertID)
amr.free.all()
