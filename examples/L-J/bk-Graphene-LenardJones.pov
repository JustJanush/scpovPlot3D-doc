/****h* ScPovPlot3D/IsoSurf.inc
* PURPOSE
*   Text object extensions, 
*     image:./imgs/H2O-01.jpg
*     
*   Fig.[TextExt] Example of visualisation of coulombean potential around
*                  the water molecule
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
#declare FileVer=" 2012-10-27";
#include ".\Plot3Dv3\CoordsSys.inc"
#include ".\Plot3Dv3\Potential.inc"

// actual scene definition
// #declare Photons=true;
global_settings{
  ambient_light color rgb 1 // sets an overall brightness/ambient light level in the scene
  max_trace_level 155                 // sets the maximum ray tracing bounce depth (1 or more) [5]
  assumed_gamma 2.2
  #if (Photons)                     // global photon block try to tune to Your needs
    photons {
      spacing 0.02               // specify the density of photons
      media 100
      count 100000               // alternatively use a total number of photons
      jitter 1.0                 // jitter phor photon rays
      max_trace_level 460         // optional separate max_trace_level
      adc_bailout 1/255          // see global adc_bailout
      autostop 0                 // photon autostop option
      radius 10                  // manually specified search radius
      expand_thresholds 0.2, 40  // (---Adaptive Search Radius---)
    }
  #end        
 
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
}

//#declare fac = 1.0;//E-4;
//#declare kqcoef = 1000*qeff*fac; // [A] zawiera mno¿nik 

// local in this file                         
//#declare BoxMin = 10*min_extent(WaterMolecule) + 1E-4*<-2,-2,-2>;
//#declare BoxMax = 10*max_extent(WaterMolecule);
//#declare BoxMax = <BoxMax.x,0, BoxMax.z> + 1E-4*< 2, 0, 2>;                           

// import molecules
#declare epsiqeff  = 10000;
#declare HOdist    = 1.0;
#declare RHyd      = HOdist/5;  
#declare RmLJ      = 2*RHyd;   

//================[ My new structure ]===============
#declare MaxN      = 10;
#declare MaxB      = 13;
#declare PTube = array[MaxN+1][8]  // charge, x,y,z, R, R_m, epsilon, multifunction see below
//                             [..][7] color component: [..][0]>0 Red; <0 Blue; (0;1)-hue component of HSL color; >=1-element number in Mendelejew table, may encode radius as well, 1-H, 2=He, ....
                 { // POVRay coordinates         
                   //          1      2     3     4       5     6     7       8
                           { MaxN,       0, 0,       0,    0,    0,        0, 0}, // MaxN=[0][0] - number of molecules in the table:[1..MaxN]
                           { -1, -2*HOdist, 0,  HOdist, RHyd, RmLJ, epsiqeff, 0}, //  1                       6 - 7 - 8 - 9 - 10 /+/ |  
                           { -1,   -HOdist, 0,  HOdist, RHyd, RmLJ, epsiqeff, 0}, //  2  q, x,y,z             |   |   |   |    |     |
                           { -1,         0, 0,  HOdist, RHyd, RmLJ, epsiqeff, 0}, //  3  q, x,y,z             1 - 2 - 3 - 4 -  5 /-/ |
                           { -1,    HOdist, 0,  HOdist, RHyd, RmLJ, epsiqeff, 0}, //  4                                              V +Z
                           { -1,  2*HOdist, 0,  HOdist, RHyd, RmLJ, epsiqeff, 0}, //  5
                         //====================================// q, x,y,z
                           {  1, -2*HOdist, 0, -HOdist, RHyd, RmLJ, epsiqeff, 0}, //  6
                           {  1,   -HOdist, 0, -HOdist, RHyd, RmLJ, epsiqeff, 0}, //  7  q, x,y,z
                           {  1,         0, 0, -HOdist, RHyd, RmLJ, epsiqeff, 0}, //  8  q, x,y,z
                           {  1,    HOdist, 0, -HOdist, RHyd, RmLJ, epsiqeff, 0}, //  9
                           {  1,  2*HOdist, 0, -HOdist, RHyd, RmLJ, epsiqeff, 0}, // 10
                 };        

#declare Bonds = array[MaxB+1][3] 
                 { 
                           {MaxB, 0, 0}, // holds number of bonds [1..MaxB]
                           {1,    6, 0}, // nums of two atoms, then bondage type (single, double, ionic, covalent, ...)
                           {2,    7, 0},
                           {3,    8, 0},
                           {4,    9, 0},
                           {5,   10, 0},
                         //=============\\
                           {1,    2, 2},
                           {2,    3, 2},
                           {3,    4, 2},
                           {4,    5, 2},
                         //=============\\ 
                           {6,    7, 2},
                           {7,    8, 2},
                           {8,    9, 2},
                           {9,   10, 2},
                 };


SetMoleculeZoxom(.2)                                                         
#declare Molecule = object{CreateStruct( Graphene_Atoms, Graphene_Bonds, 0.1)}; // PTube, Bonds, .05)};
object{Molecule} 

//====[ Control box ]==========
#declare dy = 0.0;  // BoxMax modifier, [0.0..>1] designates upper cross section plane posiotion: 0=XY plane, 1=no change >1 increase BoxMax
#declare BoxMin     = min_extent(Molecule);   
#declare BoxMax     = max_extent(Molecule);
#declare _MoleculeCenter = (BoxMin+BoxMax)/2; // fresh BoxMin/Max

#declare BoxMin     = _MoleculeCenter+ 2.5*(BoxMin-_MoleculeCenter);   
#declare BoxMax     = _MoleculeCenter+ 2.5*(BoxMax-_MoleculeCenter);
                                       
#declare BoxMin   =  <BoxMin.x, BoxMin.y*6, BoxMin.z>;
//                                       ^----- tune this coeeficient manually!
#declare BoxMax   =  <BoxMax.x, dy*BoxMax.y, BoxMax.z>;
//                               ^----- ..as well as this one

#if (Dodatki) // add switch: "Declare=Dodatki=1" to command line
  sphere{ BoxMin, RHyd texture{ pigment{color 2*Blue} finish{emission 1} }}
  sphere{ BoxMax, RHyd texture{ pigment{color 2*Red } finish{emission 1} }}
  box{ BoxMin, BoxMax 
       pigment{color rgbt<.75, .75, 0.75, .3>} 
       finish{reflection metallic diffuse .3 ambient .3 phong .3}
}           
#end

//=============================
#declare VVV  = CreateVLJ(Graphene_Atoms)  //PTube)
#declare VVV3 = CreateVLJ2(Graphene_Atoms) //PTube)

//#declare VVV  = CreateVFC(Graphene_Atoms)   // VVV() - without scaling factor
//#declare VVV3 = CreateVFC2(Graphene_Atoms)  // VVV() - WITH scaling factor

SetAccuracy(.01)
SetMaxGrad(5)

    SetIsoTexture( texture{ pigment{color rgbt<.429, .429, .429, .13>}} )
    #declare ISO3 = MakeEquiPlane( VVV, 0, BoxMin, BoxMax ) // zero plane
    object{ISO3}
    
//    SetIsoTexture( texture{ pigment{color rgbt<0,.1,1,.13>}} )
//    #declare ISO1 = MakeEquiPlane( VVV, -0.035/500000, BoxMin, BoxMax )
//    //object{ISO1}
    
//    SetIsoTexture( texture{ pigment{color <0,.1,1>}} )
//    #declare ISO2 = MakeEquiPlane( VVV, -0.07/500000, BoxMin, BoxMax )
//    //object{ISO2}
    
    SetIsoTexture( texture{ pigment{color rgbt<1,0.1,0,.13>}} )
    #declare ISO4 = MakeEquiPlane( VVV, 125/25, BoxMin, BoxMax )
    object{ISO4}
    
    SetIsoTexture( texture{ pigment{color rgbt<1,0.1,.5, .13>}} )
    #declare ISO5 = MakeEquiPlane( VVV, 10000/25, BoxMin, BoxMax )
    object{ISO5}
    
    SetIsoTexture( texture{ pigment{color rgbt<1, 1, 0.1,.15>}} )
    #declare ISO6 = MakeEquiPlane( VVV, -1.0/25, BoxMin, BoxMax )
    object{ISO6}

    SetIsoTexture( texture{ pigment{color rgbt<1, .5, 0.1,.15>}} )
    #declare ISO6 = MakeEquiPlane( VVV, -2.0/25, BoxMin, BoxMax )
    object{ISO6}

    SetIsoTexture( texture{ pigment{color rgbt<.1, .1, 1,.15>}} )
    #declare ISO7 = MakeEquiPlane( VVV, -1/50, BoxMin, BoxMax )
    object{ISO7}
    
    //PotentialMap(  VVV3, y, .01,     20 )      
    //PotentialLines(VVV3, y, .0, 111, 5 )
    //PotentialLines(VVV3, x, .0, 20, 500000 )
    //PotentialLines(VVV3, x, -1, 20, 500000 )
    //PotentialLines(VVV3, x, -2, 20, 500000 )
    //PotentialLines(VVV3, x, -4, 20, 500000 )
//======================================                              
// InsertCartesianArrows_LD(5, .03, .4)   
background {rgbt <.21, .21, .21, 1>}

light_source { < -1,  30, -1> 1 }
light_source { <  1,  30, -1> 1 }

light_source { < -1,  30,  1> 1 }
light_source { <  1,  30,  1> 1 }

light_source { <   0, -10,   0> 1 }


//#declare BoxMin     = min_extent(ISOx);   
//#declare BoxMax     = max_extent(ISOx);
//#declare _MoleculeCenter = (BoxMin+BoxMax)/2; // fresh BoxMin/Max
                                                                        
                                                                        
SetCameraTarget( _MoleculeCenter.x, _MoleculeCenter.y, _MoleculeCenter.z-.2)           
IntelligentEyeT( 13, 40, 8)       
