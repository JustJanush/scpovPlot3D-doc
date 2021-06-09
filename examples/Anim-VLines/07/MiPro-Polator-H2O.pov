/****h* ScPovPlot3D/Interpolator.inc
* PURPOSE
*   Interpolation of scalar field by means of linear approximation, see Lit: [?]
*   consolidates varius vol. render techniques for MiPro Conference.
*     image:./imgs/H2O-01.jpg
*     
*   Fig.[Interpolator] Example of interpolation inside bunch of cubes 4now
*   |html <hr width=50% align="left">
*     *********************************************************
*     **   Tested on PovRay 3.7.                             **
*     **   License: GNU GPL                                  **
*     **   Homepage:    http://scpovplot3d.sourceforge.net   **
*     *********************************************************
*     **   version: 3.1.RC1 (& have a nice time ;)             **
*     *********************************************************
* AUTHOR
*  Janusz Opi³a Ph.D.
*   jmo@agh.edu.pl, janusz.opila@gmail.com
*   Dept. of Applied Informatics
*   AGH University of Science & Technology, Cracow, Poland
*   Maintained by Janusz Opi³a Ph.D.
* COPYRIGHT
*  GNU GPL v.3 License
*  (c) 2012-now by Janusz Opi³a Ph.D.         
*  AGH University of Science and Technology
*
****
  end of RoboDoc comment
****/

#version 3.7;
#declare FileVer=" 2015-02-05";                
// Set Library Path to +LC:\Users\mirfak\Documents\POV-Ray\Sceny\dev02\examples
// +Q11  +UA +W800 +H800  +RP5  +FN8 +A +AM2 +R3   Declare=Variant=1 Declare=Debug=0 Declare=radio=0 Declare=Photons=0  +LC:\Users\mirfak\Documents\POV-Ray\Sceny\dev02\examples\Plot3Dv3
#include "colors.inc"  // import CHSL2RGB macro
#include ".\Plot3Dv3\CommonDefs.inc"
#include ".\Plot3Dv3\Cameras.inc"

#include ".\Plot3Dv3\scFinish.inc"
#include ".\Plot3Dv3\CoordsSys.inc"
#include ".\Plot3Dv3\Potential.inc"     

// actual scene definition
// #declare Photons=true;
global_settings{
  ambient_light color rgb 1      // sets an overall brightness/ambient light level in the scene
  max_trace_level 155            // sets the maximum ray tracing bounce depth (1 or more) [5]
  assumed_gamma 2.2
  #if (Photons)                  // global photon block try to tune to Your needs
  #debug "Photons switched - ON. Relax and drink a coffee."
    photons {
      spacing 0.02               // specify the density of photons
      media 100
      count 100000               // alternatively use a total number of photons
      jitter 1.0                 // jitter phor photon rays
      max_trace_level 460        // optional separate max_trace_level
      adc_bailout 1/255          // see global adc_bailout
      autostop 0                 // photon autostop option
      radius 10                  // manually specified search radius
      expand_thresholds 0.2, 40  // (---Adaptive Search Radius---)
    }
  #end        
 
 #if (radio) 
 #debug "Radiosity switched - ON. Relax and take a shower..."
  radiosity{
      pretrace_start 0.08
      pretrace_end   0.01
      count 150
      nearest_count 10
      error_bound 0.5
      recursion_limit 3
      low_error_factor 0.5
      gray_threshold 0.0
      minimum_reuse 0.005
      maximum_reuse 0.2
      brightness 1
      adc_bailout 0.005
  }         
  #end        
}
#declare dd  = 0.00025;        
#declare dd3 = pow(dd,3);

#declare _Accura = 0.001;      
#macro SetAccura(_Acr) // setter
  #declare _Accura = _Acr; 
#end 

            
#ifndef (Variant) // 4Safety
    #declare Variant  = 0;
#end


#macro DrawQVBox(_QVr)  // default look
  merge{    // bottom face
    sphere{ <_QVr[0][0][0][2], _QVr[0][0][0][3], -_QVr[0][0][0][1]>, 2*rr 
        texture{ pigment{color rgb <0, 0, 0.8>}} finish{scDullMirror ambient .7}}
    
    sphere{ <_QVr[N1-1][0][0][2], _QVr[N1-1][0][0][3], -_QVr[N1-1][0][0][1]>, rr 
        texture{ pigment{color rgb 0.25}} finish{scDullMirror ambient .7}}
    
    sphere{ <_QVr[0][N2-1][0][2], _QVr[0][N2-1][0][3], -_QVr[0][N2-1][0][1]>, rr 
        texture{ pigment{color rgb 0.25}} finish{scDullMirror ambient .7}}
    
    sphere{ <_QVr[N1-1][N2-1][0][2], _QVr[N1-1][N2-1][0][3], -_QVr[N1-1][N2-1][0][1]>, rr 
        texture{ pigment{color rgb 0.25}} finish{scDullMirror ambient .7}}
    
    // upper face
    sphere{ <_QVr[0][0][N3-1][2], _QVr[0][0][N3-1][3], -_QVr[0][0][N3-1][1]>, rr 
        texture{ pigment{color rgb 0.55}} finish{scDullMirror ambient .7}}
    
    sphere{ <_QVr[N1-1][0][N3-1][2], _QVr[N1-1][0][N3-1][3], -_QVr[N1-1][0][N3-1][1]>, rr 
        texture{ pigment{color rgb 0.55}} finish{scDullMirror ambient .7}}
    
    sphere{ <_QVr[0][N2-1][N3-1][2], _QVr[0][N2-1][N3-1][3], -_QVr[0][N2-1][N3-1][1]>, rr 
        texture{ pigment{color rgb 0.55}} finish{scDullMirror ambient .7}}
    
    sphere{ <_QVr[N1-1][N2-1][N3-1][2], _QVr[N1-1][N2-1][N3-1][3], -_QVr[N1-1][N2-1][N3-1][1]>, 2*rr 
        texture{ pigment{color rgb <0.8,0,0>}} finish{scDullMirror ambient .7}}
    
    // X-X axis
    cylinder{<_QVr[   0][0][0][2], _QVr[   0][0][0][3], -_QVr[   0][0][0][1]>, 
             <_QVr[N1-1][0][0][2], _QVr[N1-1][0][0][3], -_QVr[N1-1][0][0][1]> rc
                texture{ pigment{color rgb <.8,.2,.2>}} finish{scDullMirror ambient .7}
    }
    
    cylinder{<_QVr[   0][N2-1][0][2], _QVr[   0][N2-1][0][3], -_QVr[   0][N2-1][0][1]>, 
             <_QVr[N1-1][N2-1][0][2], _QVr[N1-1][N2-1][0][3], -_QVr[N1-1][N2-1][0][1]> rc
                texture{ pigment{color rgb <.8,.2,.2>}} finish{scDullMirror ambient .7}
    }
    
    cylinder{<_QVr[   0][0][N3-1][2], _QVr[   0][0][N3-1][3], -_QVr[   0][0][N3-1][1]>, 
             <_QVr[N1-1][0][N3-1][2], _QVr[N1-1][0][N3-1][3], -_QVr[N1-1][0][N3-1][1]> rc
                texture{ pigment{color rgb <.8,.2,.2>}} finish{scDullMirror ambient .7}
    }
    
    cylinder{<_QVr[   0][N2-1][N3-1][2], _QVr[   0][N2-1][N3-1][3], -_QVr[   0][N2-1][N3-1][1]>, 
             <_QVr[N1-1][N2-1][N3-1][2], _QVr[N1-1][N2-1][N3-1][3], -_QVr[N1-1][N2-1][N3-1][1]> rc
                texture{ pigment{color rgb <.8,.2,.2>}} finish{scDullMirror ambient .7}
    }
    
    
    
    // Y-Y axis
    cylinder{<_QVr[0][0][   0][2], _QVr[0][0][   0][3], -_QVr[0][0][   0][1]>, 
             <_QVr[0][0][N3-1][2], _QVr[0][0][N3-1][3], -_QVr[0][0][N3-1][1]> rc
                texture{ pigment{color rgb <.2,.2,.9>}} finish{scDullMirror ambient .7}
    }
    
    cylinder{<_QVr[N1-1][N2-1][   0][2], _QVr[N1-1][N2-1][   0][3], -_QVr[N1-1][N2-1][   0][1]>, 
             <_QVr[N1-1][N2-1][N3-1][2], _QVr[N1-1][N2-1][N3-1][3], -_QVr[N1-1][N2-1][N3-1][1]> rc
                texture{ pigment{color rgb <.2,.2,.9>}} finish{scDullMirror ambient .7}
    }
    
    cylinder{<_QVr[N1-1][0][   0][2], _QVr[N1-1][0][   0][3], -_QVr[N1-1][0][   0][1]>, 
             <_QVr[N1-1][0][N3-1][2], _QVr[N1-1][0][N3-1][3], -_QVr[N1-1][0][N3-1][1]> rc
                texture{ pigment{color rgb <.2,.2,.9>}} finish{scDullMirror ambient .7}
    }
    
    cylinder{<_QVr[0][N2-1][0][2], _QVr[0][N2-1][   0][3], -_QVr[0][N2-1][   0][1]>, 
             <_QVr[0][N2-1][0][2], _QVr[0][N2-1][N3-1][3], -_QVr[0][N2-1][N3-1][1]> rc
                texture{ pigment{color rgb <.2,.2,.9>}} finish{scDullMirror ambient .7}
    }
             
             
    // xYP - axis         
    cylinder{<_QVr[0][0][0][2], _QVr[0][   0][0][3], -_QVr[0][   0][0][1]>, 
             <_QVr[0][N2-1][0][2], _QVr[0][N2-1][0][3], -_QVr[0][N2-1][0][1]> rc
                  texture{ pigment{color rgb <.2,.9,.2>}}
                  finish{scDullMirror ambient .7}
    }
    
    cylinder{<_QVr[N1-1][0][0][2],    _QVr[N1-1][   0][0][3], -_QVr[N1-1][   0][0][1]>, 
             <_QVr[N1-1][N2-1][0][2], _QVr[N1-1][N2-1][0][3], -_QVr[N1-1][N2-1][0][1]> rc
                  texture{ pigment{color rgb <.2,.9,.2>}}
                  finish{scDullMirror ambient .7}
    }
    
    cylinder{<_QVr[0][   0][N3-1][2], _QVr[0][   0][N3-1][3], -_QVr[0][   0][N3-1][1]>, 
             <_QVr[0][N2-1][N3-1][2], _QVr[0][N2-1][N3-1][3], -_QVr[0][N2-1][N3-1][1]> rc
                  texture{ pigment{color rgb <.2,.9,.2>}}
                  finish{scDullMirror ambient .7}
    }
    cylinder{<_QVr[N1-1][   0][N3-1][2], _QVr[N1-1][   0][N3-1][3], -_QVr[N1-1][   0][N3-1][1]>, 
             <_QVr[N1-1][N2-1][N3-1][2], _QVr[N1-1][N2-1][N3-1][3], -_QVr[N1-1][N2-1][N3-1][1]> rc
                  texture{ pigment{color rgb <.2,.9,.2>}}
                  finish{scDullMirror ambient .7}
    }
  } // merge  
#end         


  
   //==========================================================//
  //                                                          //
 // show me the structure - visualization of the whole grid  //
//==========================================================//
#macro DrawAllCells(_QVr) // default look
    merge{                             
        // All QV nodes
        #for(in1,0,N1-1)   
           #for(in2,0,N2-1)
              #for(in3,0,N3-1)                
                 sphere{ <_QVr[in1][in2][in3][2],_QVr[in1][in2][in3][3], -_QVr[in1][in2][in3][1]>, rr  
                    texture{ pigment{color rgb 0.8*<in1/N1, in2/N2, in3/N3>+0.2}}
                    finish{scDullMirror ambient .7}
                 }
              #end
           #end
        #end                  
        
        // cylinders along xR - axis - RED
        union{
            #for(inz,0,N3-1)
               #for(iny,0,N2-1)
                        cylinder{<_QVr[0]   [iny][inz][2], _QVr[0]   [iny][inz][3], -_QVr[0]   [iny][inz][1]>, 
                                 <_QVr[N1-1][iny][inz][2], _QVr[N1-1][iny][inz][3], -_QVr[N1-1][iny][inz][1]> rc
                        texture{ pigment{color rgb <.8,.2,.2>}}
                        finish{scDullMirror ambient .7}
                     }
               #end
            #end  
        }                
        
        
        // cylinders along yRr - axis - GREEN
        union{
            #for(inz,0,N3-1)
               #for(inx,0,N1-1)
                        cylinder{<_QVr[inx][   0][inz][2], _QVr[inx][   0][inz][3], -_QVr[inx][   0][inz][1]>, 
                                 <_QVr[inx][N2-1][inz][2], _QVr[inx][N2-1][inz][3], -_QVr[inx][N2-1][inz][1]> rc
                        texture{ pigment{color rgb <.2,.9,.2>}}
                        finish{scDullMirror ambient .7}
                     }
               #end
            #end  
        }                
        
        // cylinders along zRr - axis - BLUE
        union{
            #for(iny,0,N2-1)
               #for(inx,0,N1-1)
                        cylinder{<_QVr[inx][iny][   0][2], _QVr[inx][iny][   0][3], -_QVr[inx][iny][   0][1]>, 
                                 <_QVr[inx][iny][N3-1][2], _QVr[inx][iny][N3-1][3], -_QVr[inx][iny][N3-1][1]> rc
                        texture{ pigment{color rgb <.2,.2,.9>}}
                        finish{scDullMirror ambient .7}
                     }
               #end
            #end  
        }                
    } // merge
#end

///*   
   //====================================//
  //  key interpolation function        //
 //  internal POVRay coordinates!!!    // 
//====================================//
#declare VV1 = function(sx, sy, sz,   // xP, yP, zP
                        zx, zy, zz,   // lower, left, front corner of elementary cell
                        Q1d, Q2d, Q3d, Q4d,  // loads on boottom nodes
                        Q1u, Q2u, Q3u, Q4u,  // loads on upper nodes
                        ddd, ddd3){ 
// new version
  (
      (// bottom hex face  
         (-Q1d*(-sz+zz-ddd)-Q3d*(-zz+sz)) *(zx+ddd-sx) //Q13d*(1-z) +..reversed z-axis!
         +
         (-Q2d*(-sz+zz-ddd)-Q4d*(-zz+sz)) *(sx-zx)    //Q24d*z 
      )*(zy+ddd-sy)
      +
      (// upper hex face
         (-Q1u*(-sz+zz-ddd)-Q3u*(-zz+sz)) *(zx+ddd-sx) //Q13u*(1-z) +..reversed z-axis!
         +
         (-Q2u*(-sz+zz-ddd)-Q4u*(-zz+sz)) *(sx-zx)    //Q24u*z 
      )*(sy-zy)
  )/ddd3

} 
//*/

   //===========================//
  //                           //
 //    Drawing facility       //
//===========================//
#macro DrawInterpolSrf(_trsh, _dd) // _trsh - isosurface threshold
#debug concat("=>dd= ",str(_dd,6,4),"==\n")   
#declare ttr_ = concat(" Trsh= ", str(_trsh,6,4));   

#declare _dd3 = pow(_dd,3);
    #for(izn,0,N3-2)
        #debug concat(" Processing layer.. [", str(izn,3,0) "] by treshold", ttr_ ," \n")    
        #for(iyn,0,N2-2)
            #for(ixn,0,N1-2)                               
                // do not analyze "empty" cells
                #if (_trsh<min(
                                    QVr[  ixn][  iyn][izn][0],   // Q1d
                                    QVr[  ixn][iyn+1][izn][0],   // Q2d
                                    QVr[ixn+1][  iyn][izn][0],   // Q3d
                                    QVr[ixn+1][iyn+1][izn][0],   // Q4d
                                    
                                    // upper face
                                    QVr[  ixn][  iyn][izn+1][0], // Q1u
                                    QVr[  ixn][iyn+1][izn+1][0], // Q2u
                                    QVr[ixn+1][  iyn][izn+1][0], // Q3u
                                    QVr[ixn+1][iyn+1][izn+1][0]  // Q4u
                           ) 
                    | _trsh>max(
                                    QVr[  ixn][  iyn][izn][0],   // Q1d
                                    QVr[  ixn][iyn+1][izn][0],   // Q2d
                                    QVr[ixn+1][  iyn][izn][0],   // Q3d
                                    QVr[ixn+1][iyn+1][izn][0],   // Q4d
                                    
                                    // upper face
                                    QVr[  ixn][  iyn][izn+1][0], // Q1u
                                    QVr[  ixn][iyn+1][izn+1][0], // Q2u 
                                    QVr[ixn+1][  iyn][izn+1][0], // Q3u
                                    QVr[ixn+1][iyn+1][izn+1][0]  // Q4u
                            )      
                     |  (QVr[ixn][iyn][izn][4])     
                    ) // do nothing
//                #debug concat(ttr_ " Cell [", str(ixn,3,0), "Cell [", str(ixn,3,0),"][", str(iyn,3,0),"][",str(izn,3,0),"] ",str(QVr[ixn][iyn][izn][4],3,0), "  ommited, \n")    
                #else
                    isosurface{
                      function{VV1( x,y,z,    // conversion from real to internal representation
                                    QVr[ixn][iyn][izn][2], QVr[ixn][iyn][izn][3], -QVr[ixn][iyn][izn][1],      // "root" of elementary cell
                                    
                                    // bottom face
                                    QVr[  ixn][  iyn][izn][0],   // Q1d
                                    QVr[  ixn][iyn+1][izn][0],   // Q2d
                                    QVr[ixn+1][  iyn][izn][0],   // Q3d
                                    QVr[ixn+1][iyn+1][izn][0],   // Q4d
                                    
                                    // upper face
                                    QVr[  ixn][  iyn][izn+1][0], 
                                    QVr[  ixn][iyn+1][izn+1][0],  
                                    QVr[ixn+1][  iyn][izn+1][0], 
                                    QVr[ixn+1][iyn+1][izn+1][0],
                                    _dd, _dd3
                                    ) 
                      }               
                      contained_by{ 
                            box{
                                 <QVr[  ixn][  iyn][  izn][2], QVr[  ixn][  iyn][  izn][3], -QVr[  ixn][  iyn][  izn][1]>,     // "root" of elementary cell >,
                                 <QVr[ixn+1][iyn+1][izn+1][2], QVr[ixn+1][iyn+1][izn+1][3], -QVr[ixn+1][iyn+1][izn+1][1]>      // "root" of elementary cell
                            }
                      }
                      threshold _trsh
                      accuracy _Accura 
                      all_intersections
                      max_gradient 1.1  
                      evaluate .1, 1.4, .7
                      open               
                      texture{_IsoTexture}
                    }
                #end
            #end 
        #end 
    #end
#end


//===============[ End of INC file ]==============\\


SetMoleculeZoom(.6) // ***                                                        

//#declare Molecule = object{CreateStruct( Graphene_Atoms, Graphene_Bonds, 0.2)}; // PTube, Bonds, .05)};
 #declare Molecule = object{CreateStruct( H2O_Atoms, H2O_Bonds, 0.07)}; // *** PTube, Bonds, .05)};7
 object{Molecule} 

//====[ Control box ]==========
#declare dy = .950;  // *** BoxMax modifier, [0.0..>1] designates upper cross section plane position: 0=XY plane, 1=no change >1 increase BoxMax
#declare BoxMin     = min_extent(Molecule);   
#declare BoxMax     = max_extent(Molecule);
#declare _MoleculeCenter = (BoxMin+BoxMax)/2; // fresh BoxMin/Max

#ifdef (AddOns)
    #if (AddOns) // add switch: "Declare=Dodatki=1" to command line
      sphere{ BoxMin, RHyd texture{ pigment{color 2*Blue} finish{emission 1} }}
      sphere{ BoxMax, RHyd texture{ pigment{color 2*Red } finish{emission 1} }}
      box{ BoxMin, BoxMax 
           pigment{color rgbt<.75, .75, 0.75, .3>} 
           finish{reflection metallic diffuse .3 ambient .3 phong .3}
    }           
    #end         
#else
    #warning "Add switch: \"Declare=AddOns=1\" to command line to see private stuff ;) \n"       
#end

#declare BoxMin = _MoleculeCenter+ 3.4*(BoxMin-_MoleculeCenter);   
#declare BoxMax = _MoleculeCenter+ 3.8*(BoxMax-_MoleculeCenter);
//                                 ^----- tune these coeeficients manually!
                                       
#declare BoxMin =  <BoxMin.x, 3*BoxMin.y, BoxMin.z>; // *** 
//                            ^----- tune this coeeficient manually!
#declare BoxMax =  <.0+1*BoxMax.x, (1-dy)*BoxMax.y, BoxMax.z>; // ***
//                                  ^----- ..as well as this one

SetAccuracy(.005)//.1)     
SetMaxGrad(1)

// ***** BEGIN
#declare VVV_LJ1 = CreateVLJ(Graphene_Atoms)     // *** PTube)
#declare VVV_LJ2 = CreateVLJ2(Graphene_Atoms, 2) // *** PTube)
#declare VVV_C1  = CreateVFC(H2O_Atoms)          // *** VVV() - without scaling factor
#declare VVV_C2  = CreateVFC2(H2O_Atoms, 1E9)    // *** VVV() - WITH scaling factor

/* This does not work!!
#declare FreeFun = function(x,y,z) {
        #local c=x;
        #local d=z;
         d
} 
*/

#declare VVV = VVV_C2; // *** wa¿ne!!
// ***** END


// available "Variant" values
#declare _IsoSurf  = 1;
#declare _Interpol = 2;
#declare _TextreFunc = 3;
#declare _PseudoPart = 4;

#declare RGB1 = color_map { 
          [ 0.05 rgbt < 1, 1, 1, .0>]
          [ 0.95 rgbt < 0, 0, 0, .0>]
          [ 0.95 rgbt < 1, 1, 1, .0>]
          };     


#declare RGB2 = color_map { 
          [ 0.00 rgbt < .2, .2, .2, .0>]
          [ 0.05 rgbt < .5, .5, .5, .0>]
          [ 0.15 rgbt < 0, 0, 1, .0>]
          [ 0.45 rgbt < 0, .7, 0, .0>]
          [ 0.75 rgbt < 1, .7, 0, .0>]
          [ 0.95 rgbt < 1, 0, 0, .0>]
          };     

// *** NEW!
#declare RGB3 = color_map { 
          [ 0.000 rgbt < 1 , 1 , 1 , .0>]
          [ 0.010 rgbt < .4, .4, .4, .0>]
          [ 0.050 rgbt < .4, .4, .6, .0>]
          [ 0.200 rgbt < 0, 0, 1, .0>]
          [ 0.300 rgbt < 0, .5, 1, .0>]
          [ 0.400 rgbt < 0, 1, .50, .0>]
          [ 0.500 rgbt < 0, 1, 0, .0>]
          [ 0.600 rgbt < 1, 1, 0, .0>]
          [ 0.700 rgbt < 1, 1, 0, .0>]
          [ 0.800 rgbt < 1, 0, 0, .0>]
          [ 0.900 rgbt < 1, 0, 0.7, .0>]
          [ 0.990 rgbt < 1, 0, 1, .0>]
          [ 0.990 rgbt < 2, 2, 2, .0>]
          [ 0.999 rgbt < 2, 2, 2, .0>]
          [ 0.999 rgbt < 1, 0, 1, .0>]
          [ 1.000 rgbt < 1, 0, 1, .0>]
          };     
          
#declare _ColorMap = RGB3;

cylinder{ <.0+0.3*BoxMax.x, 0.0001, .8*BoxMax.z>, 
          <.0+0.3*BoxMax.x, .9999, .8*BoxMax.z>,
          .08 // *** // ***, 1, .2
    texture{ 
              pigment{ 
                gradient y
                color_map{_ColorMap}
              }                               
                finish{ scDullMirror }
    }     
    scale y*3
//    translate 4*x
}



#switch (Variant)
////===================
    #case ( _IsoSurf )
        // tune BoxMin/Max to this vis style
//SetIsoTexture( texture{ pigment{color rgbt<.5, .5, .5, .0>}} )

    #local gMin = 0;//-3.5;
    //#local gMin2= .64;
    //#local gMax = 2.05;                                         
    #local gMin2 = 0;
    #local gMax = 1900;
    
    SetGradientAccuracy(.0001)      

    #local grV0 = fn_Gradient(VVV);     // *** VVV
    // #local grV = function{ clip(log(grV0(x,y,z)-gMin)-gMin2, 0, gMax-gMin2)/(gMax-gMin2) };
    #local grV = function{ clip(grV0(x,y,z)-gMin2, 0, gMax-gMin2)/(gMax-gMin2) };
      SetIsoTexture( texture{ 
          pigment{ 
            
            // *** fn_Gradient(function{ clip(log(VVV(x,y,z)-gMin),0,gMax)/(gMax-gMin) }) 
            // *** function{ clip(log(VVV(x,y,z)), 0, gMax )/gMax})}
            // *** fn_Gradient(function{ clip((VVV(x,y,z)-gMin),0,gMax-gMin)/(gMax-gMin) }) 
            // *** function{ clip(log(VVV(x,y,z)), 0, gMax )/gMax})}
            function{ grV(x,y,z) }  //function{ clip(log(VVV(x,y,z)), 0, gMax )/gMax})}

            color_map{_ColorMap}
          }
          
          finish{ Dull } // *** scDullMirror }
      } )
          
          
// SetIsoTexture( texture{ pigment{color rgbt<.5,0.5,0.5,.0>}} )
             #declare ISO3 = MakeEquiPlane( VVV, 0, BoxMin, BoxMax ) // *** zero plane
       object{ISO3}
        
//      SetIsoTexture( texture{ pigment{color rgbt<1,0.1,0,.0>}} )
       #declare ISO4 = MakeEquiPlane( VVV, 25.0, BoxMin, BoxMax )   // ***
       object{ISO4}

//      SetIsoTexture( texture{ pigment{color rgbt<1,0.1,.5, .0>}} )
       #declare ISO5 = MakeEquiPlane( VVV, 50, BoxMin, BoxMax )
       object{ISO5}
    
///    SetIsoTexture( texture{ pigment{color rgbt<1, .5, 0.1,.15>}} )
    #declare ISO6 = MakeEquiPlane( VVV, -25.0, BoxMin, BoxMax )
    object{ISO6}

///    SetIsoTexture( texture{ pigment{color rgbt<.1, .1, 1,.15>}} )
    #declare ISO7 = MakeEquiPlane( VVV, -50.0, BoxMin, BoxMax )
    object{ISO7}

//      SetIsoTexture( texture{ pigment{color rgbt<1, 1, 0.0,.0>}} )
//       #declare ISO6 = MakeEquiPlane( VVV, -1000.95, BoxMin, BoxMax ) 
       //object{ISO6}
    
      #break 
////===================

    #case (_Interpol)  
// compute node values
// option: read them in

        #declare N1  = 150;
        #declare N2  = 100;
        #declare N3  = 30;              
        #declare dd  = .1;        
        #declare dd3 = pow(dd,3);
        #declare rr  = .3*dd;
        #declare rc  = rr/3;
        #declare QVr = array[N1][N2][N3][5];
        // QV stores [xR index][yR index][zR index][V,x,y,z]real!!! \
        // real data convert to internal representation: x>-zp, y>xp, z>yp
        // no translation                                        
        // bottom   1  2  3
        // the view" from top:
        // bottom layer
        // zR +--> yR     Q1--Q2
        //    |           |   |
        //    V xR        Q3--Q4
        // upper layer  
        // zR +--> yR     Q5--Q6
        //    |           |   |
        //    V xR        Q7--Q8
        //============================//
        
        // Definition of the grid - usually will be completed outside POVRay
        // then read in by another #macro
        // format of the file: 
        // "Data_title"
        // "Data_Comment,
        // "X_Label", "Y_Label", "Z_Label",
        // node-node interval - 'dd', N(rows in file), Nx, Ny, Nz (array dimensions),
        // ix, iy, iz, Q, flag,
        // .
        // .                               
        // .
        
        #for(in1,0,N1-1)                                                 
           #for(in2,0,N2-1)
              #for(in3,0,N3-1)              
                    // real coordinates!!  
                    #local xxx = in1*dd-BoxMax.z-1.5; 
                    #local yyy = in2*dd+BoxMin.x-1.5; 
                    #local zzz = in3*dd+BoxMax.y-3.0;
                    #local rrV = xxx*xxx+yyy*yyy+zzz*zzz;
        
                    #declare QVr[in1][in2][in3][0] = VVV(yyy, zzz, -xxx); //zzz*(xxx-yyy)/(1+rrV); // (xxx*xxx-yyy*yyy-zzz*zzz)/(.1+rrV); //exp(-(rrV*rrV)/2); 
                    #declare QVr[in1][in2][in3][1] = xxx; 
                    #declare QVr[in1][in2][in3][2] = yyy;
                    #declare QVr[in1][in2][in3][3] = zzz;
                    #declare QVr[in1][in2][in3][4] = 0*mod(in1, 10); // 0-render, !=0 do not render or another action (possibly will be rendered with another resolution)
        
              #end
           #end
        #end                  
        
// volume marker - prism edges        
        DrawQVBox(QVr)                 
        #declare rr  = .1*dd;
        #declare rc  = rr/3;
        // alternatively draw ALL cells (a lot of mess!)
        // DrawAllCells(QVr)                  
                                           
      SetIsoTexture( texture{ pigment{color rgbt<.5, .5, .5, 0>}} )
        DrawInterpolSrf( .0, dd)

      SetIsoTexture( texture{ pigment{color rgbt<1,0.1,0,0>}} )
        DrawInterpolSrf( 500, dd)

      SetIsoTexture( texture{ pigment{color rgbt<1,0.1,.5, 0>}} )
        DrawInterpolSrf( 20000, dd)
        
      SetIsoTexture( texture{ pigment{color rgbt<1, 1, 0., .0>}} )
        DrawInterpolSrf( -.95, dd)
        
#ifdef (AddOns)
    #if (AddOns) // add switch: "Declare=Dodatki=1" to command line
        // marks center of the "Molecule"
        sphere{<(N2-2)*dd/2, 0*(N3-1)*dd/2, -(N1-1)*dd/2> .5 texture{pigment{color rgb <1,0,0>}} }    
    #end         
#else
    #warning "Add switch: \"Declare=AddOns=1\" to command line to see private stuff ;) \n"       
#end
 
        #break                       

    #case (_TextreFunc)  
       #warning "Vis. by texture map. mostly implemented :) :/ \n" // 4safety        
    #local gMin = -3;
    #local gMax = 1e5;
    SetGradientAccuracy(.0001)      

    SetRGBMap(RGB2)
    
    //PotentialMapG(  fn_Gradient(function{clip(log(VVV(x,y,z)), 0, gMax )/gMax}), y, .0 )      
    // PotentialMapG(  function{clip(log(VVV(x,y,z)-gMin), gMin, gMax )/(gMax-gMin)}, y, .0 )      
     PotentialMapG(function{ clip( VVV(x,y,z)-gMin, 0, gMax-gMin)/(gMax-gMin)}, y, .0 )      
                                
    // PotentialMap(  function{log(VVV(x,y,z))}, x, .0 )          
    // PotentialLines(function{log(VVV(x,y,z)+3)}, y, .0, 1 ) // inline anonymous function
    // PotentialLines(function{log(VVV(x,y,z)+3)}, x, .0, 1 )
    // PotentialLines(function{-log(VVV(x,y,z)+3)}, y, 4.5-clock*9, 2 )
//    PotentialLines(function{log(VVV(x,y,z)+3)}, x, 9.6-clock*12.6, 2 )

    #break        
    #case (_PseudoPart) 
        #warning "Vis. using pseudo-particles. Not fully implemented, sorry :/ \n" // 4safety       
        #local _num =10000;  
        #local _shot=0.0;    
        #local _vmin=1e11;
        #local _vmax=0.0;
        #declare ax = 12; // x range
        #declare ay = 3;  // y range
        #declare az = 16; // sd(z) 0h55m01sat 3.5 - nice ;)              
        
                              
        #declare R1 = seed(1);                                            
        #declare R2 = seed(69);                      
        #declare R3 = seed(val(datetime(now)));                      
    
        #local Vtrsh=0.0001;                  
        #local rr=1;
    
    #local grd = 0.0;    // ***
    #local maxgrd = grd; // ***
    
    #local gMin = 0;
    #local gMax = 010.0600;
    SetGradientAccuracy(.0001)      

        #while (_num>0)  
           #local xr = ax*rand(R1)-3; // povray coords
           #local yr = ay*rand(R2)-3;  
           #local zr = az*rand(R3)-8;
           
           #local Vval= VVV(xr, yr, zr); // ***       
           //#if (mod(_num,100)=0)         // ***
           //  #debug concat( "| num= ", str(_num,7,0), "| rx= ",str(-zr,10,3), "| rx= ",str(xr,10,3), "| rx= ",str(yr,10,3),  " | V(r)= ",str(Vval,12,9) "\n") // ***
           //#end                          // ***
           
           #if (abs(Vtrsh-Vval)<rr)                         
               #local p0 = <xr,yr,zr>;              // cell found
               // sphere{p0, .04 pigment{color rgb 1 }} //test sphere

                                          
               #local vGrd   = vGradient(function{VVV(x,y,z)}, p0);  // compute gradient in the center of the cell 'p0'
               #local Lvgr   = vlength(vGrd); // VALUE of the grad   
               #local _vmin  = min(Lvgr, _vmin);// choose NEW "min"
               #local _vmax  = max(Lvgr, _vmax);// choose NEW "max"
               #local nvGrd  = vnormalize(vGrd);    // ..normalize it, then ..
               #if(Lvgr>01.01)
                  #local p0 = p0+(Vtrsh-VVV(p0.x, p0.y, p0.z))*nvGrd/Lvgr; // ..compute and apply correction to the point 'p0'
               #end
               
               // #local grd   = Gradient_Length(function{ (clip(VVV(x,y,z)-gMin, 0, gMax-gMin)/(gMax-gMin)) }, p0); 
               #local grd    = clip(Gradient_Length(VVV, p0)-gMin, 0, gMax-gMin)/(gMax-gMin);
               #local maxgrd = max(maxgrd,grd);         // ***  
               #local pgmnt  = pigment{ function{grd} color_map{RGB2}};
               sphere{
                    p0, clip(abs(Vtrsh-VVV(p0.x, p0.y, p0.z)), .02, .04) //max(0.05, rr/50)
                    texture{ pgmnt }
               }
               #local _num=_num-1;
          #end                    
          #local _shot=_shot+1;
          #if (mod(_shot,10000)=0)
             #debug concat("n= ",str(_shot,9,0), "| num= ", str(_num,7,0), " | Vmin= ",str(_vmin,10,3), "| Vgrad= ",str(grd,10,3), " | Vmax= ",str(_vmax,10,3),  "\n") // ***
          #end
        #end                                  
        #debug concat("\n=========\n| tot.= ",str(_shot,10,0), " |\n=========")
        #debug concat("\n| Vmin= ",str(_vmin,10,3), "| Max_grd= ",str(maxgrd,10,5),"| Vmax= ",str(_vmax,10,3)," |\n=========") // ***
    #break
    #else 
      #warning "No variant chosed, use Declare=Variant=# on the command line\n"        
      #warning "where '#' stands for 1, 2, 3 or 4 by now:\n"        
      #warning " 1 - raw isosurface \n"        
      #warning " 2 - raw interpolation \n"
      #warning " 3 - texture \n"
      #warning " 4 - pseudo particles \n"
#end  // switch(Variant

// Lights, background & Camera:

// ***    background {rgbt 1*<1, 1, 1, 1>}
    background {rgbt <0, 0, 0, 1>}

// /* 
    light_source { < 0, 10,  0> 1 shadowless}

    light_source { < -10,   7, -1> 1 shadowless}
    light_source { <  10,   7, -1> 1 shadowless}
    light_source { < -10, -10,  3> 1 shadowless}
    light_source { <  10, -10,  3> 1 shadowless}
// */ 
/* 
    light_source { < 0,  10,  0> 1 shadowless}
    light_source { < 0, -10,  0> 1 shadowless}

    light_source { < -5,   7, -3> 1 shadowless}
    light_source { <  3,   7, -3> 1 shadowless}
    light_source { < -5,   7,  1> 1 shadowless}
    light_source { <  3,   7,  1> 1 shadowless}
*/
//   InsertCartesianArrows_LDWT(1.2, .02, 1.2, 1.2,
//        texture{pigment{color rgb 1} finish { scDullMirror } }) // visualization ofcoords system                 

//    InsertCartesianArrows_LD(9, .15, .4)   // ***

    SetCameraTarget( _MoleculeCenter.x, _MoleculeCenter.y-1.20, _MoleculeCenter.z+.5)           
    IntelligentEyeT( 85, 55, 11)       
//CameraLights(15, 4, .1, 60, 30 )
    
