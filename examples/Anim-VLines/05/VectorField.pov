/****E* ScPovPlot3D/VectorFields.pov
* PURPOSE                     
*   This module is in introductory state.
*   It contains macros for representation of vector field in form of set of vectors symbolised by
*   various shapes and color coding systems. Direction and strength of vector at given space point
*   can be visualised in different ways. While direction of the vector can be shown by main axis
*   of some figure, for example cone or cylinder, strength and turn can be represented by length 
*   or color or volume or so on. It depends mainly on the goal of visualisation. 
*   I think, that representation of vectors in single plane is most informative by now.
*   Besides that  we need superimposition of source objects, as coils, charges, permanent magnet 
*   poles or even oceanic bed if one takes into account visualisation of oceanic currents.
*     image:./Plot3dv3/imgs/VectorField.png
*     |html <hr width=50% align="left">
*     *********************************************************
*     **   Tested on PovRay 3.7                              **
*     **   License: GNU GPL                                  **
*     **   Homepage:    http://scpovplot3d.sourceforge.net   **
*     *********************************************************
*     **   version: 3.0.6 (& have a nice time ;)             **
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
#include "colors.inc"
#include "textures.inc"                  
#include "Plot3Dv3\Cameras.inc"
#include "Plot3Dv3\ColorMaps.inc"
#include "Plot3Dv3\Potential.inc"


#version 3.7;
#local _FileName   = "VectorFields.pov";
#local _FileVerMaj = "3.1.0.1";
#local _FileVerMin = "2015-12-29";
#debug concat("\n[==> ", _FileName,", ver: ", _FileVerMaj, ", build:  ", _FileVerMin, " <==]\n")
      
// some defaults             
#local NP = 1;      
#declare NPmax= NP*NP*NP;        
                           
#macro DeclareVectors(_VS, _NP)
    #declare NPmax= _NP*_NP*_NP;        
    #declare SpcPntsVecs = array[NPmax][3]; // [][0] - space point, [][1] field vector direction, [][2] field strength
    #for(iz, 0, _NP-1)
        #for(iy, 0, _NP-1)
            #for(ix, 0, _NP-1)                                                  
                #local indx = _NP*_NP*iz+_NP*iy+ix ;
                #declare SpcPntsVecs[indx][0] = <ix+.2, iy+.2, iz+.2>;
                #declare VL = vlength(SpcPntsVecs[indx][0]);                                    
                #declare SpcPntsVecs[indx][1] = SpcPntsVecs[indx][0]/VL;
                #declare SpcPntsVecs[indx][2] = <_VS(QQ, VL),0,0>; // THIS Line needs refinment           
            #end   
        #end                      
    #end   
#end

/*
"Field measured on Hawaii on 2012-12-12",
111, 
x, y, z, Vx, Vy, Vz // position <x,y,z>_REAL! and vector strength <vx, vy, vz>_REAL 
x, y, z, Vx, Vy, Vz // real!! ......
*/
#macro ImportVectors(_filenam)

   #fopen DataFile _filenam read
   #read(DataFile, _Descript)
   #read(DataFile, _NumPts)


        #for(iz, 0, NP-1)
            #for(iy, 0, NP-1)
                #for(ix, 0, NP-1)                                                  
                    #local indx = NP*NP*iz+NP*iy+ix ;
                    #declare SpcPntsVecs[indx][0] = <ix+.2, iy+.2, iz+.2>;
                    #declare VL = vlength(SpcPntsVecs[indx][0]);                                    
                    #declare SpcPntsVecs[indx][1] = SpcPntsVecs[indx][0]/VL;
                    #declare SpcPntsVecs[indx][2] = <_VS(QQ, VL),0,0>; // THIS Line needs refinment           
                #end   
            #end                      
        #end   

    #fclose

#end

#macro DrawVectorField()                        
    #declare VF = union{
        #for(SpcPnt, 0, NPmax-1)
           #declare VL = vlength(SpcPntsVecs[SpcPnt][1]);                                                 
           cone{<0, 0, 0>, Rr, SpcPntsVecs[SpcPnt][1], 0 
                  texture{ 
                     pigment {
                         gradient SpcPntsVecs[SpcPnt][1]
                         color_map {
                              [0.00 color rgbt<.0, .0, .13, .0>]
                              [0.30 color rgbt<.0, .0, .33, .0>]
                              [0.33 color rgbt<.7, .7, .70, .0>]
                              [0.36 color rgbt<.7, .0, .00, .0>]
                              [1.00 color rgbt<.7, .0, .00, .0>]
                         }        
                         scale 1.05*VL                        
                         translate -0.025*SpcPntsVecs[SpcPnt][1]/VL 
                      }
                      finish{ phong .3 reflection .5 diffuse .3 ambient .9 }
                  }                               
                  translate -0.33*SpcPntsVecs[SpcPnt][1]/VL               
                  scale ScFac*SpcPntsVecs[SpcPnt][2].x                                
                  translate SpcPntsVecs[SpcPnt][0]                  
           }
    
        #end   
    };                                            

    object{VF}
#end                                 

global_settings{
  ambient_light color rgb <.9,.9,.9> // sets an overall brightness/ambient light level in the scene
  max_trace_level 5                 // sets the maximum ray tracing bounce depth (1 or more) [5]
  assumed_gamma 1.0  
  
  #if (radio)
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

  #if (Photons)                  // global photon block try to tune to Your needs
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
}

// from old H2O file

SetMoleculeZoom(.6) // ***                                                        

//#declare Molecule = object{CreateStruct( Graphene_Atoms, Graphene_Bonds, 0.2)}; // PTube, Bonds, .05)};
 #declare Molecule = object{CreateStruct( H2O_Atoms, H2O_Bonds, 0.07) }; // *** PTube, Bonds, .05)};7
 object{Molecule } 

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
#declare VVV_LJ1 = CreateVLJ(Graphene_Atoms)    // *** PTube)
#declare VVV_LJ2 = CreateVLJ2(Graphene_Atoms, 2) // *** PTube)
#declare VVV_C1  = CreateVFC(H2O_Atoms)          // ***  VVV() - without scaling factor
#declare VVV_C2  = CreateVFC2(H2O_Atoms, 1E9)       // *** VVV() - WITH scaling factor

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
          [ 0.000 rgbt < 1.0, 1.0, 1.0, 0.0>]
          [ 0.010 rgbt < 0.4, 0.4, 0.4, 0.0>]
          [ 0.050 rgbt < 0.4, 0.4, 0.6, 0.0>]
          [ 0.200 rgbt < 0.0, 0.0, 1.0, 0.0>]
          [ 0.300 rgbt < 0.0, 0.5, 1.0, 0.0>]
          [ 0.400 rgbt < 0.0, 1.0, 0.5, 0.0>]
          [ 0.500 rgbt < 0.0, 1.0, 0.0, 0.0>]
          [ 0.600 rgbt < 1.0, 1.0, 0.0, 0.0>]
          [ 0.700 rgbt < 1.0, 1.0, 0.0, 0.0>]
          [ 0.800 rgbt < 1.0, 0.0, 0.0, 0.0>]
          [ 0.900 rgbt < 1.0, 0.0, 0.7, 0.0>]
          [ 0.990 rgbt < 1.0, 0.0, 1.0, 0.0>]
          [ 1.000 rgbt < 1.0, 1.0, 0.0, 0.0>]
          };     
          

//=======================

light_source { < 0, 10,     0> 1 shadowless}
light_source { < -10,   7, -1> 1 shadowless}
light_source { <  10,   7, -1> 1 shadowless}
light_source { < -10, -10,  3> 1 shadowless}
light_source { <  10, -10,  3> 1 shadowless}

background {rgb <.10, .10, .10>/7}

#declare RR = .25;  // promieñ punktu g³ównego
#declare Rr = .20;  // promieñ punktu kontrolnego
#declare rr = .08;  // promieñ walca ³¹cz¹cego punkty g³ówne i kontrolne
#declare CardinalPoint = texture{pigment {color    Red} finish{ambient .36 reflection metallic phong .3 emission .2}}
#declare ControlPoint  = texture{pigment {color Orange} finish{ambient .36 reflection metallic phong .3 emission .2}}
#declare ControlLine   = texture{pigment {color Silver} finish{ambient .6 reflection metallic phong .3}}
#declare Trans =.55; 

#local QQ    = 1;
#local ScFac = .75;                                                                             

// some userspace functions 
#declare VSlog = function(_Q, _L){log(_Q/_L+1)}
#declare VSn   = function(_Q, _L){1}

//DeclareVectors(VSn, 15)
//DrawVectorField()
//=================[ Epilog ]====================

//#declare TP =  (min_extent(VF)+max_extent(VF))/2;                                                          

#local gMin = 1; // 0;//-3.5;
//#local gMin2= .64;
//#local gMax = 2.05;                                         
#local gMin2 = 1.; //0
#local gMax = 3.5;  //1900

#local grV0 = fn_Gradient(VVV);     // *** VVV
 #local grV = function{ clip(log(grV0(x,y,z)-gMin)-gMin2, 0, gMax-gMin2)/(gMax-gMin2) };
//#local grV = function{ clip(grV0(x,y,z)-gMin2, 0, gMax-gMin2)/(gMax-gMin2) };
#declare _ColorMap = RGB3;

SetIsoTexture( 
     texture{ 
          pigment{ 
            
            // *** fn_Gradient(function{ clip(log(VVV(x,y,z)-gMin),0,gMax)/(gMax-gMin) }) 
            // *** function{ clip(log(VVV(x,y,z)), 0, gMax )/gMax})}
            // *** fn_Gradient(function{ clip((VVV(x,y,z)-gMin),0,gMax-gMin)/(gMax-gMin) }) 
            // *** function{ clip(log(VVV(x,y,z)), 0, gMax )/gMax})}
            function{ grV(x,y,z) }  
            color_map{_ColorMap}
          }
          finish{ Dull } // *** scDullMirror }
    } 
)

/****M* VectorField.inc/VectorCF
* PURPOSE
*  this macro computes Coulomb vector field from 
*  assembly of charged molecules, passed as global variable _Qs
*  
* SYNOPSIS
*/
#macro VectorCF ( _rr, _sc )
/* 
* INPUTS
*  3Dvector _rr  - 3D space point, units - [A] 
*  float    _sc  - additional scaling factor, may be left _sc==1
*  garray   _PTs - array[N+1][8] of charges, passed as global variable, use macro "SetChargeTable(H2O_Atoms)", refer Potential.inc for details
* SEE ALSO
*  None
*   
******/

// _PTs( charge, x,y,z, [... not used ...] )    
 #local _sc = _sc*qeff;
 #local _vx = 0.0;
 #local _vy = 0.0;
 #local _vz = 0.0;

 #for(i, 1, _PTs[0][0]) 
     #local _R3 = pow(pow(_rr.x-_PTs[i][1],2) + pow(_rr.y-_PTs[i][2],2) + pow(_rr.z-_PTs[i][3],2), 1.5); // R^3
     #local _vx = _vx + _PTs[i][0]*(_rr.x-_PTs[i][1])/_R3;                     
     #local _vy = _vy + _PTs[i][0]*(_rr.y-_PTs[i][2])/_R3;
     #local _vz = _vz + _PTs[i][0]*(_rr.z-_PTs[i][3])/_R3;
 #end                       
 #local _vx = _vx*_sc;
 #local _vy = _vy*_sc;
 #local _vz = _vz*_sc;
 <_vx, _vy, _vz>
#end  // end  VectorCF() macro

#declare Graphite = rgb 0.10;

#declare dx=1.0/3*<1, 0, 0> ;
#declare dy=<0,0,0>; // float
#declare dz=1.0/3.0*<0, 0,-1> ;
#declare rr  =.025;
#declare rrc =.035;
#declare rrr = 20000;

#local Dalp  = 0; // degrees
#local Dalp  = radians(Dalp); // radians
#local radii = 1.2*HOdist;
#local ymin   = 0;
#local ymax   = .5;  
#local yy    = ymin+clock*(ymax-ymin);

#macro DeclareStructure(_Atoms)
  #declare _PTs = _Atoms;
#end                            

DeclareStructure(H2O_Atoms)                                                     

#local Nlin = 6;                 
#local Npts = 160;
#local StrmLineTable = array[Nlin+1][Npts+1][4]; // [0][0][0] - num of streamlines, [0][*][*] - reserved
                                                 // [k][0][0] - num of points in `k`-th streamline, may depend on field's structure n<=Npts   

#local StrmLineTable[0][0][0] = Nlin;
#local _h = .0329; // Runge-Kutty step
// control points

#local alpha0 = -5.0; // degrees
#local alpha0 = radians(alpha0); // radians
#declare _P0 = radii*<cos(alpha0),        yy,  -sin(alpha0)>;  // test: cylinder{P00,VV,StrmLineTable[1][0]*1000 texture{pigment{color rgb<1.0, .3,  .3>} finish{phong .7 reflection metallic}}}

// lets make streamlines
#local StrmLineTable[1][0][0] = Npts;
#local _VV = VectorCF(_P0, 1); // local electric vector [V/m]
#local StrmLineTable[1][1][0] = rr; //rrr*vlength(VV); // vector length [V/m]=>radius of given node of streamline
#local StrmLineTable[1][1][1] = _P0.x; // [A]
#local StrmLineTable[1][1][2] = _P0.y; // [A]
#local StrmLineTable[1][1][3] = _P0.z; // [A]

#for(ind, 2, Npts) // Runge-Kutta of 4 order           
    #local _vn = vnormalize( _VV ); // unary vector parallel to local vector field
    #local _k1 = _h*_vn;            // k1 (R-K)
    #local _VV = VectorCF((<_P0.x, _P0.y, _P0.z>+_h/2.0*_vn), 1); // x+h/2
    #local _k2 = _h*vnormalize( _VV );
    #local _k3 = _k2;
    #local _VV = VectorCF((<_P0.x, _P0.y, _P0.z>+_h*_vn), 1); // x+h/2
    #local _k4 = _h*vnormalize( _VV ); 

    #local _P1 = _P0+(_k1+2*_k2+2*_k3+_k4)/6.0;
     
    #local StrmLineTable[1][ind][0] = rr; //rrr*vlength(VV); // vector length [V/m]=>radius of given node of streamline
    #local StrmLineTable[1][ind][1] = _P1.x; // [A]
    #local StrmLineTable[1][ind][2] = _P1.y; // [A]
    #local StrmLineTable[1][ind][3] = _P1.z; // [A]
    #local _P0 = _P1;
    #local _VV = VectorCF( _P0, 1); // local electric vector [V/m]

#end


// 2nd line                            
#local alpha0 = 109.45; // degrees
#local alpha0 = radians(alpha0); // rad
#declare _P0 = radii*<cos(alpha0),        yy,  -sin(alpha0)>;  // test: cylinder{P00,VV,StrmLineTable[1][0]*1000 texture{pigment{color rgb<1.0, .3,  .3>} finish{phong .7 reflection metallic}}}

// lets make streamlines
#local StrmLineTable[2][0][0] = Npts;
#local _VV = VectorCF(_P0, 1); // local electric vector [V/m]
#local StrmLineTable[2][1][0] = rr; //rrr*vlength(VV); // vector length [V/m]=>radius of given node of streamline
#local StrmLineTable[2][1][1] = _P0.x; // [A]
#local StrmLineTable[2][1][2] = _P0.y; // [A]
#local StrmLineTable[2][1][3] = _P0.z; // [A]

#for(ind, 2, Npts) // Runge-Kutta of 4 order           
    #local _vn = vnormalize( _VV ); // unary vector parallel to local vector field
    #local _k1 = _h*_vn;            // k1 (R-K)
    #local _VV = VectorCF((<_P0.x, _P0.y, _P0.z>+_h/2.0*_vn), 1); // x+h/2
    #local _k2 = _h*vnormalize( _VV );
    #local _k3 = _k2;
    #local _VV = VectorCF((<_P0.x, _P0.y, _P0.z>+_h*_vn), 1); // x+h/2
    #local _k4 = _h*vnormalize( _VV ); 

    #local _P1 = _P0+(_k1+2*_k2+2*_k3+_k4)/6.0;
     
    #local StrmLineTable[2][ind][0] = rr; //rrr*vlength(VV); // vector length [V/m]=>radius of given node of streamline
    #local StrmLineTable[2][ind][1] = _P1.x; // [A]
    #local StrmLineTable[2][ind][2] = _P1.y; // [A]
    #local StrmLineTable[2][ind][3] = _P1.z; // [A]
    #local _P0 = _P1;
    #local _VV = VectorCF( _P0, 1); // local electric vector [V/m]

#end

// 3rd line                            
#local alpha0 = 119.45; // degrees
#local alpha0 = radians(alpha0); // rad
#declare _P0 = radii*<cos(alpha0),        yy,  -sin(alpha0)>;  // test: cylinder{P00,VV,StrmLineTable[1][0]*1000 texture{pigment{color rgb<1.0, .3,  .3>} finish{phong .7 reflection metallic}}}

// lets make streamlines
#local StrmLineTable[3][0][0] = Npts;
#local _VV = VectorCF(_P0, 1); // local electric vector [V/m]
#local StrmLineTable[3][1][0] = rr; //rrr*vlength(VV); // vector length [V/m]=>radius of given node of streamline
#local StrmLineTable[3][1][1] = _P0.x; // [A]
#local StrmLineTable[3][1][2] = _P0.y; // [A]
#local StrmLineTable[3][1][3] = _P0.z; // [A]

#for(ind, 2, Npts) // Runge-Kutta of 4 order           
    #local _vn = vnormalize( _VV ); // unary vector parallel to local vector field
    #local _k1 = _h*_vn;            // k1 (R-K)
    #local _VV = VectorCF((<_P0.x, _P0.y, _P0.z>+_h/2.0*_vn), 1); // x+h/2
    #local _k2 = _h*vnormalize( _VV );
    #local _k3 = _k2;
    #local _VV = VectorCF((<_P0.x, _P0.y, _P0.z>+_h*_vn), 1); // x+h/2
    #local _k4 = _h*vnormalize( _VV ); 

    #local _P1 = _P0+(_k1+2*_k2+2*_k3+_k4)/6.0;
     
    #local StrmLineTable[3][ind][0] = rr; //rrr*vlength(VV); // vector length [V/m]=>radius of given node of streamline
    #local StrmLineTable[3][ind][1] = _P1.x; // [A]
    #local StrmLineTable[3][ind][2] = _P1.y; // [A]
    #local StrmLineTable[3][ind][3] = _P1.z; // [A]
    #local _P0 = _P1;
    #local _VV = VectorCF( _P0, 1); // local electric vector [V/m]

#end

// 4th line                            
#local alpha0 = -15; // degrees
#local alpha0 = radians(alpha0); // rad
#declare _P0 = radii*<cos(alpha0),        yy,  -sin(alpha0)>;  // test: cylinder{P00,VV,StrmLineTable[1][0]*1000 texture{pigment{color rgb<1.0, .3,  .3>} finish{phong .7 reflection metallic}}}

// lets make streamlines
#local StrmLineTable[4][0][0] = Npts;
#local _VV = VectorCF(_P0, 1); // local electric vector [V/m]
#local StrmLineTable[4][1][0] = rr; //rrr*vlength(VV); // vector length [V/m]=>radius of given node of streamline
#local StrmLineTable[4][1][1] = _P0.x; // [A]
#local StrmLineTable[4][1][2] = _P0.y; // [A]
#local StrmLineTable[4][1][3] = _P0.z; // [A]

#for(ind, 2, Npts) // Runge-Kutta of 4 order           
    #local _vn = vnormalize( _VV ); // unary vector parallel to local vector field
    #local _k1 = _h*_vn;            // k1 (R-K)
    #local _VV = VectorCF((<_P0.x, _P0.y, _P0.z>+_h/2.0*_vn), 1); // x+h/2
    #local _k2 = _h*vnormalize( _VV );
    #local _k3 = _k2;
    #local _VV = VectorCF((<_P0.x, _P0.y, _P0.z>+_h*_vn), 1); // x+h/2
    #local _k4 = _h*vnormalize( _VV ); 

    #local _P1 = _P0+(_k1+2*_k2+2*_k3+_k4)/6.0;
     
    #local StrmLineTable[4][ind][0] = rr; //rrr*vlength(VV); // vector length [V/m]=>radius of given node of streamline
    #local StrmLineTable[4][ind][1] = _P1.x; // [A]
    #local StrmLineTable[4][ind][2] = _P1.y; // [A]
    #local StrmLineTable[4][ind][3] = _P1.z; // [A]
    #local _P0 = _P1;
    #local _VV = VectorCF( _P0, 1); // local electric vector [V/m]

#end

// 5th line                            
#local alpha0 = 114.45; // degrees
#local alpha0 = radians(alpha0); // rad
#declare _P0 = radii*<cos(alpha0),        yy,  -sin(alpha0)>;  // test: cylinder{P00,VV,StrmLineTable[1][0]*1000 texture{pigment{color rgb<1.0, .3,  .3>} finish{phong .7 reflection metallic}}}

// lets make streamlines
#local StrmLineTable[5][0][0] = Npts;
#local _VV = VectorCF(_P0, 1); // local electric vector [V/m]
#local StrmLineTable[5][1][0] = rr; //rrr*vlength(VV); // vector length [V/m]=>radius of given node of streamline
#local StrmLineTable[5][1][1] = _P0.x; // [A]
#local StrmLineTable[5][1][2] = _P0.y; // [A]
#local StrmLineTable[5][1][3] = _P0.z; // [A]

#for(ind, 2, Npts) // Runge-Kutta of 4 order           
    #local _vn = vnormalize( _VV ); // unary vector parallel to local vector field
    #local _k1 = _h*_vn;            // k1 (R-K)
    #local _VV = VectorCF((<_P0.x, _P0.y, _P0.z>+_h/2.0*_vn), 1); // x+h/2
    #local _k2 = _h*vnormalize( _VV );
    #local _k3 = _k2;
    #local _VV = VectorCF((<_P0.x, _P0.y, _P0.z>+_h*_vn), 1); // x+h/2
    #local _k4 = _h*vnormalize( _VV ); 

    #local _P1 = _P0+(_k1+2*_k2+2*_k3+_k4)/6.0;
     
    #local StrmLineTable[5][ind][0] = rr; //rrr*vlength(VV); // vector length [V/m]=>radius of given node of streamline
    #local StrmLineTable[5][ind][1] = _P1.x; // [A]
    #local StrmLineTable[5][ind][2] = _P1.y; // [A]
    #local StrmLineTable[5][ind][3] = _P1.z; // [A]
    #local _P0 = _P1;
    #local _VV = VectorCF( _P0, 1); // local electric vector [V/m]

#end

// 6th line                            
#local alpha0 = -10; // degrees
#local alpha0 = radians(alpha0); // rad
#declare _P0 = radii*<cos(alpha0),        yy,  -sin(alpha0)>;  // test: cylinder{P00,VV,StrmLineTable[1][0]*1000 texture{pigment{color rgb<1.0, .3,  .3>} finish{phong .7 reflection metallic}}}

// lets make streamlines
#local StrmLineTable[6][0][0] = Npts;
#local _VV = VectorCF(_P0, 1); // local electric vector [V/m]
#local StrmLineTable[6][1][0] = rr; //rrr*vlength(VV); // vector length [V/m]=>radius of given node of streamline
#local StrmLineTable[6][1][1] = _P0.x; // [A]
#local StrmLineTable[6][1][2] = _P0.y; // [A]
#local StrmLineTable[6][1][3] = _P0.z; // [A]

#for(ind, 2, Npts) // Runge-Kutta of 4 order           
    #local _vn = vnormalize( _VV ); // unary vector parallel to local vector field
    #local _k1 = _h*_vn;            // k1 (R-K)
    #local _VV = VectorCF((<_P0.x, _P0.y, _P0.z>+_h/2.0*_vn), 1); // x+h/2
    #local _k2 = _h*vnormalize( _VV );
    #local _k3 = _k2;
    #local _VV = VectorCF((<_P0.x, _P0.y, _P0.z>+_h*_vn), 1); // x+h/2
    #local _k4 = _h*vnormalize( _VV ); 

    #local _P1 = _P0+(_k1+2*_k2+2*_k3+_k4)/6.0;
     
    #local StrmLineTable[6][ind][0] = rr; //rrr*vlength(VV); // vector length [V/m]=>radius of given node of streamline
    #local StrmLineTable[6][ind][1] = _P1.x; // [A]
    #local StrmLineTable[6][ind][2] = _P1.y; // [A]
    #local StrmLineTable[6][ind][3] = _P1.z; // [A]
    #local _P0 = _P1;
    #local _VV = VectorCF( _P0, 1); // local electric vector [V/m]

#end


/****M* VectorField.inc/DrawStreamLines
* PURPOSE
*  Draws streamlines based on Streamline table _SL, spline type _splt and texture _txt
* SYNOPSIS
*/
#macro DrawStreamLines ( _SL, _splt, _txt )// - table, spline, texture
/* 
* INPUTS
*  garray     _SL   - streamlines table ov N lines over N_i nodes
*  integer    _splt - interpolation type: 2 - bezier spline, 3 - cubic spline, other - linear_spline
*  texturedef _txt  - texture attached to ALL lines, may be simple or complex, ex. functional
* SEE ALSO
*  None
*   
******/
    #for(jnd, 1, _SL[0][0][0]) 
        sphere_sweep{//
            #switch(_splt)
              #case( 2 ) 
                 b_spline      
              #break
              #case( 3 ) 
                 cubic_spline  
              #break
              #else 
                 linear_spline
            #end  
           _SL[jnd][0][0]
                        
           #local _np = _SL[1][0][0];             
           #for(ind, 1, _np)
              < _SL[jnd][ind][1], _SL[jnd][ind][2], _SL[jnd][ind][3]>, _SL[1][ind][0]  
           #end              
    
           texture{ _txt   
                finish { scDullMirror } 
           } // end of texture 
        }
        sphere{ <_SL[jnd][ 1 ][1], _SL[jnd][ 1 ][2], _SL[jnd][ 1 ][3]>, rrc texture{pigment{color rgb< 1.0, .1,  .1>} finish{scDullMirror}}}                                                                  
        sphere{ <_SL[jnd][_np][1], _SL[jnd][_np][2], _SL[jnd][_np][3]>, rrc texture{pigment{color rgb<  .1, .1, 1.0>} finish{scDullMirror}}}
    #end
    
#end

//=================================
DrawStreamLines( StrmLineTable, 1, _IsoTexture )   

#local alf_ = 51;
#local alfr =radians(alf_);
cylinder{ <.0, 0.0001, .0>, 
          <.0,  .9999, .0>,
          .04 // *** // ***, 1, .2
    texture{ 
              pigment{ 
                gradient y
                color_map{_ColorMap}
              }                               
              finish{ scDullMirror emission 1 }
    }     
    scale y*2        
    translate <0,-1,0>
    rotate 90*x
    rotate (90-alf_)*y
    translate 1.5*<cos(alfr),0,-sin(alfr)>
    translate 0.5*<sin(alfr), 0,cos(alfr)>
//    translate 4*x
}

SetCameraTarget(.0, 0.0, .2)           
IntelligentEyeT( alf_, 57, 8)


