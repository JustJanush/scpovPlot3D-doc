/****E* ScPovPlot3D/MultiCrits.pov
* PURPOSE                     
*   image:./Plot3dv3/imgs/MultiCrits.png
*
*   Fig.[Shape1] Example of histogram render/1
*   |html <hr width=50% align="left"> 
*     *********************************************************
*     **   Tested on PovRay 3.7                              **
*     **   License: GNU GPL                                  **
*     **   Homepage:    http://scpovplot3d.sourceforge.net   **
*     *********************************************************
*     **   version: 3.0.7 (& have a nice time ;)             **
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
#declare _FileVerMaj ="3.0.7";
#declare _FileVerMin ="2014-12-22";
#declare _FileName   ="tst-Particles.pov";                                

#debug concat("\n[==> ", _FileName,", ver: ", _FileVerMaj, ", build:  ", _FileVerMin, " <==]\n")


// here are macros essential to Your plot 
#include "colors.inc"
#include "stones.inc"
#include ".\Plot3Dv3\scFinish" 
#include ".\Plot3Dv3\Histogram.inc" 
/*********************/
//#debug _TextFont;

// if You like long renderings times, uncomment line with "TurnPhotonsUp()":
// and use glassy textures ;)
// TurnPhotonsUp()

// should be set before scene description: this is for future extensions
// The 'clock' variable is introduced for animations, useful after adding to command line
// at least one switch  +KFF###, where ### stands for number of LAST frame. First frame is always numbered '1'
// other animation options are also available (+KC, +KFI, +KFS, +KFE, ...). Read POVRay manual.

PrepareCamera(clock*360-25, 27, 70) // (longitude[deg], latitude[deg], distance[m]): 

// These settings affects strongly behaviour of the renderer
// Feel free to adjust them to Your need
global_settings{
  ambient_light color rgb 1 // sets an overall brightness/ambient light level in the scene
  max_trace_level 25                // sets the maximum ray tracing bounce depth (1 or more) [5]
  assumed_gamma 2.20
  #if (Photons)                     // global photon block try to tune to Your needs
#warning "|==>Photons==>ON"         // Usually it is the best choice, but very slow if glass is present in the scene
    photons {
      spacing 0.02               // specify the density of photons
      media 100
      count 100000               // alternatively use a total number of photons
      jitter 1.0                 // jitter phor photon rays
      max_trace_level 450        // optional separate max_trace_level
      adc_bailout 1/255          // see global adc_bailout
      autostop 0                 // photon autostop option
      radius 10                  // manually specified search radius
      expand_thresholds 0.2, 40  // (---Adaptive Search Radius---)
    }
  #end
}

#declare txMica = texture{pigment{color Mica}};// finish{scDullMirror}}; 

/**************************************************/
// You need some light, also.
// Modify, add or remove due to your need!
light_source { <  5, 25,   -5> color rgb 1  shadowless }  // natural light source, white
light_source { < 10,  9,  -10> color rgb .51 shadowless }  // natural light source, gray
light_source { <  7,  7,  -10> color rgb .51 shadowless }  // as above
light_source { <  4,  9,  -10> color rgb .51 shadowless }  // as above
light_source { <  1,  7,  -10> color rgb .51  shadowless }  // natural light source, white
//light_source { <  15,  1   -.1> color rgb .51 shadowless  }  // natural light source, white

background {rgbt .1*<1.0, 1.0, 1.0, 0.00>}
//background {rgbt <.10, .10, .10, 1>}
/********** Histogram definitions starts here ****************/
//SetBarObject(StandardBar)
//SetRGBFTColor( .951, .10, .10, .08, .30 ) 
//SetBarObject(CylinderMQuadroBar(.315, .32, 45))
//AddBar(4.93)                         // represent vertical value 5.93
//SetBarObject(McrEllipsoConeBasic(1, 0))
//AddBar(4.22)                  
                                 
//SetBarObjectFully(EllipsoCone)        
//SetRGBFTColor(0, .33, .77, 0.08, .30)  // RGB color is set up 'Blue'
//SetBarObject(McrEllipsoConeBasic(.5, -45))
//SetBarObject(CylinderM2QuadroBar(.72, 45))
//AddBar(6.9)  

//SetInterior(inGlass1)
//SetRGBFTColor_2(  0, .53, .27, .08, .10)  // RGB color is set up 'Blue'
//SetRGBFTColor_1( 1, .9, .0, 1, 0 ) 
//SetBarObjectFullyNSc(HollowBar(3, 1.51, 30, 0.31))
//SetBarObjectFullyNSc(HollowCylinder(4, 1, 30, 0.1))
//AddBar(-2)  // silently ignores bar height as it is set in SetBarObjectFullyNSc(HollowBar(-3, 1, 30, 0.1)) macro call
//SetInterior(inGlass)
//SetBarObjectFully(EllipsoCone)
//AddBar(.23)
//SetInterior(inGlass2)                    
//SetRGBFTColor_1( .2, .90, .20, .08, .2 )    
//SetBarObjectFullyNSc(HollowBarOpen(3, 1, 30, 0.1))
//SetBarObjectFullyNSc(HollowCylinderOpen(3, 1, 30, 0.1))
//SetBarObject(CylinderBar)
//AddBar(1)        

//================================

//SetLettrSize(0.8)                         // 
//SetLettrDepth(0.2)
//SetLettrBase(-0.7)
//Set_dX_dY(1, 1)                         // offsets between bars
//Set_Xside_Yside(1.5, 1.5)                 // horizontal dimensions of the bars
//Set_ext_dX_dY(1, 1.5) 
// TurnPhotonsDown()
       
#declare mx = 1000;
#declare my = 1000;
#declare rm = .005; // radius
#declare ax = 10;  // sd
#declare ay = 10;  // sd
#declare az = 1.5; // sd               
                      

#declare R1 = seed(1);                      
#declare R2 = seed(2);                      
#declare R3 = seed(val(datetime(now)));                      

//#declare podivajtese = merge {
    #for(i,1,mx)
        #for(j,1,my)
           #local xr = ax*rand(R1); 
           #local yr = ay*rand(R2);  
           #local zr = rand(R3);  

           #if (zr<.45)
               #local clr = <.1,.1, 1-zr>;
               #local clr = vnormalize(clr);
           #elseif (zr>.55)
               #local clr = <zr,.1,0.1>;        
               #local clr = vnormalize(clr);
           #else
               #local clr = <1,1,1>;           
           #end
           
           #local zr = xr*yr/10+az*(zr-.5); 
           
           sphere{
               <yr,zr,-xr>, rm
               texture{pigment{color rgb clr}}
           }
        #end
    #end
//} 
// 
//SetYScaleFactor(.2) // Optional: vertical scale exaggeration 


//object{
//  Text("")
//  texture{ T_Copper_3C }// T_Copper_3C} //txMica}
//  scale <1.0,1.2,1>*1.3
//  translate <9.5, 7., -.30>
//}
                              
/***************************************/
// view type,
// uncomment only one of following:
SetPerspective()              // perspective projection
//SetOrthographic(16)         // ortographic projection camera is shifted by '16'  
InsertCartesianArrows_LD(11, .2, .75)

// CameraOn(2.25)                 // Prepared camera is in effect now with field of view set to 50%


//object{ podivajtese }
//#declare bmin = min_extent(podivajtese);
//#declare bmax = max_extent(podivajtese);

//SetCameraTarget((bmin.x+bmax.x)/2-.1, (bmin.y+bmax.y)/2+.3, (bmin.z+bmax.z)/2+.0)
SetCameraTarget(5, 3, -5.5)          
#debug ChangeLog
//SetCameraTarget(5, 2, -5)
IntelligentEyeT(-26, 27, 44)

#if (Photons)
  #warning "Photons==>ON"     // Just info
#else
  #warning "Photons==>OFF"    // Just info
#end                                                                                  
                                                                                       
#debug "*********************************" 
#debug "*********************************" 

// Also possible:
// SetCameraTarget(XX, YY, ZZ)  // Camera looks at <XX, YY, ZZ> point
// CameraSpher(lam, phi, r)     // Camera is located in spherical coordinates
// CameraSpher0(lam, phi, r)    // Camera is located in spherical coordinates and looks at <0,0,0>
// IntelligentEye(25, 30, 15)   // Camera tries to catch out whole scene
// IntelligentEyeT(lam, phi, r) // Camera tries to catch out whole scene, positioned against spherical coordinates    