                               // Persistence of Vision Ray Tracer Scene Description File
// File: .pov
// Vers: 3.6
// Desc:
// Date:
// Auth:
global_settings { 
  assumed_gamma 1.0 
  number_of_waves 3 
  max_trace_level 15 
  ambient_light <.5, .5, .5>}


// various glass finishes, colors and interiors
#include "glass.inc"

// Standard pre-defined colors
#include "colors.inc"
// macros for generating and manipulating text strings
// 
#include "strings.inc"

// a lot of stone textures
// T_Stone1 - T_Stone44
#include "stones.inc"

// several different gold colors, finishes and textures
#include "golds.inc"

// various metal colors, finishes and textures
// brass, copper, chrome, silver
#include "metals.inc"     
#include "textures.inc"
        

#declare I_Glass2 =                    //Use with Bead
    interior{
       fade_distance .5              // only for this scene
       fade_power 2
       ior 2.45
       caustics 1
       fade_color <1, 2, 1>
    }
 
#declare T_GlassLight = texture {
   pigment { color red 2 green 2 blue 2 filter 0.95 }
   finish {
      ambient 0.0
      diffuse 0.0
      reflection 0.1
      phong 0.2
      phong_size 100
//      metallic
   }
}
 
// perspective (default) camera
camera {
  location  <5.0, 6.0, -10.0>
  look_at   <0.0, 0.0,  0.0>
  right     x*image_width/image_height
}


// general light definition
light_source { <10, 10, -10>  color rgb 1.0  shadowless }
light_source { <15, 20, -10>  color rgb .7 }
light_source { <-5, 20, -10>  color rgb .5 }
light_source { <5, -20, 10>  color rgb .5 }  
union {
merge {// k¹t > 180

        difference {
                merge{ // pó³ciastko
                        cylinder { -1*y, 1*y, 2}
                        torus {2, 1}                                     
                        scale <2, 1, 2>   
                }           
        
                plane {1*z, 0}
                plane {1*z, 0 rotate -120*y} // k¹t wycinka       
        }     

        difference {
                merge{
                        cylinder {-1*y,  1*y,  2}
                                        
                        torus {2, 1 }
                        scale <2, 1, 2>                           
                }                
                plane {-1*z, 0}
//                plane {-1*z, 0 rotate 45*y}         // dobre dla wycinka < 180
                        
        }                

        texture {T_Stone2}       
        finish {ambient 1.5 phong .2 }       
        rotate y*30
}


merge {// k¹t < 180

        difference {
                merge{
                        cylinder {-1*y,  1*y,  2}
                        torus {2, 1 }
                        scale <2, 1, 2>                           
                }                
                plane {-1*z, 0}
                plane {-1*z, 0 rotate 90*y}         // dobre dla wycinka < 180
        }                

        texture {T_Stone6}       
        finish {ambient 1.5 phong .2 }       
        rotate y*180
}

merge {// k¹t < 180
        difference {
                merge{
                        cylinder {-1*y,  1*y,  2}
                        torus {2, 1 }
                        scale <2, 1, 2>                           
                }                
                plane {-1*z, 0}
                plane {-1*z, 0 rotate 150*y}         // dobre dla wycinka < 180
        }                
  texture { Polished_Brass }        
//  texture { T_GlassLight }
        interior { I_Glass2 }          
        finish {ambient 1.5 phong .2 }       
        rotate y*-150

}      

rotate (clock*360)*y


}//union END