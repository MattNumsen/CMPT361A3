#include <stdio.h>
#include <GL/glut.h>
#include <math.h>
#include "global.h"
#include "sphere.h"
#include <iostream>

//
// Global variables
//
extern int win_width;
extern int win_height;

extern GLfloat frame[WIN_HEIGHT][WIN_WIDTH][3];  

extern float image_width;
extern float image_height;

extern Point eye_pos;
extern float image_plane;
extern RGB_float background_clr;
extern RGB_float null_clr;

extern Spheres *scene;
extern Spheres *board_black;
extern Spheres *board_white;

// light 1 position and color
extern Point light1;
extern float light1_intensity[3];

// global ambient term
extern float global_ambient[3];

// light decay parameters
extern float decay_a;
extern float decay_b;
extern float decay_c;

extern int shadow_on;
extern int reflection_on;
extern int refraction_on;
extern int stochastic_ref_on;
extern int supersamp_on;
extern int step_max;


extern int raysShot;
extern int raysHit;

//Vector chessboard_normal = {1.0, 1.0, 0.0};
Vector chessboard_normal = {1.0, 1.0, 0.0};
Point chessboard_point = {0.0, 5.0, 6.0};
extern int plane_on;

RGB_float black = {0, 0, 0};
RGB_float white = {1, 1, 1};

/////////////////////////////////////////////////////////////////////

/*********************************************************************
 * Phong illumination - you need to implement this!
 *********************************************************************/
RGB_float phong(Point q, Vector toEye, Vector surf_norm, Spheres *sph);
Vector getReflectedRay(Vector l, Vector n)
{
  float ndotl= vec_dot(n, l);
  Vector reflected = vec_minus(vec_scale(n, 2*ndotl), l);
  return reflected;
}

Vector getRefractedRay(Point* firstHit, Vector ray, Point* secondHit, Spheres* sph){

	//get the refracted ray based on math detailed here: http://graphics.stanford.edu/courses/cs148-10-summer/docs/2006--degreve--reflection_refraction.pdf
	// See equations (22) and (23) specifically

	Vector firstNorm = sphere_normal(*firstHit, sph);
	float N1 = 1.0;
	float N2 = sph->mat_refract;
	//std::cout<<"N2:"<<N2<<" ";
	Vector incident = vec_scale(ray, -1);
	float IdotN = -vec_dot(incident, firstNorm);

	float n, c2;


	n = N1/N2;
	c2 = sqrt(1 - (pow(n, 2)*(1 - pow(IdotN, 2))));
	Vector innerRefraction = vec_plus(vec_scale(incident, n), vec_scale(firstNorm, n*IdotN - c2));
	normalize (&innerRefraction);
	//now that we have the inverted ray, we have to track it to the point where it leaves the sphere
	
	float a, b, c, cHelper;
  	Vector shiftedCenter = get_vec(sph->center, *firstHit); //(firstHit - sph.center)
  	cHelper = vec_len(shiftedCenter);

  	a = vec_dot(innerRefraction,innerRefraction);
  	b = 2*vec_dot(innerRefraction, shiftedCenter);
  	c = cHelper*cHelper - sph->radius*sph->radius;
  	float t, t1, t2;

  	float discriminant = b*b - 4*a*c;

  	t1 = (-b + sqrt(discriminant))/(2.0*a);
    t2 = (-b - sqrt(discriminant))/(2.0*a); //one of these will be zero. The other one should be positive. 

    t = std::max(t1, t2);

    secondHit->x = firstHit->x + t*innerRefraction.x;
    secondHit->y = firstHit->y + t*innerRefraction.y;
    secondHit->z = firstHit->z + t*innerRefraction.z;

    Vector secondNorm = vec_scale(sphere_normal(*secondHit, sph), -1);
    normalize(&secondNorm); //just incase?
    IdotN = -vec_dot(innerRefraction, secondNorm);
    n = N2/N1;
    c2 = sqrt(1 - (pow(n, 2)*(1 - pow(IdotN, 2))));
    Vector secondRefraction = vec_plus(vec_scale(innerRefraction, n), vec_scale(secondNorm, n*IdotN - c2));
    normalize(&secondRefraction);
    //std::cout<<"IncidentRay (x, y, z): ("<<incident.x<<", "<<incident.y<<", "<<incident.z<<")"<<" secondRefraction (x, y, z): ("<<secondRefraction.x<<", "<<secondRefraction.y<<", "<<secondRefraction.z<<")"<<std::endl;

    return secondRefraction;

}

bool intersect_plane(Point Q, Vector V, Point Pl, Vector Pn, Point* hit){

	float t, numerator, denominator;

	numerator = Pn.x*(Q.x - Pl.x) + Pn.y*(Q.y - Pl.y) + Pn.z*(Q.z - Pl.z);
	denominator = vec_dot(V, Pn);
	//std::cout<<"Vector: {"<<V.x<<","<<V.y<<","<<V.x<<"}"<<"Normal: {"<<Pn.x<<","<<Pn.y<<","<<Pn.x<<"}"<<std::endl;;

	if (denominator != 0 || numerator == 0){ //points originating on the plane 
		t = numerator/denominator;
		//std::cout<<t<<std::endl;
		if (t < 0.0 || t > 10000) {
			return false;
		}
		hit->x = Q.x + t*V.x;
		hit->y = Q.x + t*V.y;
		hit->z = Q.x + t*V.z;

	} else { // if the denominator is zero, that means that the vector and plane are ... parallel? in some way? I think?
		return false;
	}

	//now bound it

	float x_mid = 0;
	float z_mid = -7;
	float offset = 8;
	float x_right_lim = x_mid + offset;
	float x_left_lim = x_mid - offset;
	float z_back_lim = z_mid + offset;
	float z_front_lim = z_mid - offset;
	if (hit->x <= x_right_lim && hit->x >= x_left_lim && hit->z <= z_back_lim && hit->z >= z_front_lim){
		return true;
	} else {
		return false;
	}
}

RGB_float getPlaneColor(Point hit, Vector toEye) 
{

	//first way i thought to do this is as quadrants around the point in xz plane, 0,0
	
	float x_mid = 0;
	float z_mid = -7;
	float offset = 8;
	bool blck = false;
	bool whte = false;
	float x_right_lim = x_mid + offset;
	float x_left_lim = x_mid - offset;
	float z_back_lim = z_mid + offset;
	float z_front_lim = z_mid - offset;


	RGB_float color;
	if (hit.x <= x_right_lim && hit.x >= x_mid && hit.z <= z_back_lim && hit.z >= z_mid){
		    if ((int) hit.x %2 == 0){ //even
		    	if ((int) hit.z %2 == 0){
		    	 	blck = true;
		    	} else {
		    		whte = true;
		    	}
		    } else {
		    	if ((int) hit.z %2 == 0){
		    	 	whte = true;
		    	} else {
		    		blck = true;
		    	}
		    }
		} else if (hit.x >= x_left_lim && hit.x < x_mid && hit.z <= z_back_lim && hit.z >= z_mid) {
			if ((int) hit.x %2 == 0){ //even
		    	if ((int) hit.z %2 == 0){
		    	 	whte = true;
		    	} else {
		    		blck = true;
		    	}
		    } else {
		    	if ((int) hit.z %2 == 0){
		    	 	blck = true;
		    	} else {
		    		whte = true;
		    	}
		    }
		} else if (hit.x <= x_right_lim && hit.x >= x_mid && hit.z >= z_front_lim && hit.z < z_mid) {
			if ((int) hit.x %2 == 0){ //even
		    	if ((int) hit.z %2 == 0){
		    	 	blck = true;
		    	} else {
		    		whte = true;
		    	}
		    } else {
		    	if ((int) hit.z %2 == 0){
		    	 	whte = true;
		    	} else {
		    		blck = true;
		    	}
		    }
		} else if (hit.x >= x_left_lim && hit.x < x_mid && hit.z >= z_front_lim && hit.z < z_mid) {
			if ((int) hit.x %2 == 0){ //even
		    	if ((int) hit.z %2 == 0){
		    	 	whte = true;
		    	} else {
		    		blck = true;
		    	}
		    } else {
		    	if ((int) hit.z %2 == 0){
		    	 	blck = true;
		    	} else {
		    		whte = true;
		    	}
		    }
		}


		if (blck == true){
			color = phong(hit, toEye, chessboard_normal, board_black);
		} else {
			color = phong(hit, toEye, chessboard_normal, board_white);
		}

		Vector toLight = get_vec(hit, light1);
		if (shadow_on && detectShadow(hit, toLight, scene)){
			color = clr_scale(color, 0.5);
		}
	return color;
}

RGB_float phong(Point q, Vector toEye, Vector surf_norm, Spheres *sph) {
//
// do your thing here
//

  //need a few things first. We already have v and n, but need vec to the light source (l) and the reflected ray vector (r)
  Vector toLight = get_vec(q, light1); 
  float delta = vec_len(toLight); //delta is used in the attenuation coefficient calculation
  normalize(&toLight);


  Vector reflected = getReflectedRay(toLight, surf_norm);
  normalize(&reflected);

  float attenuation = 1/(decay_a + decay_b*delta + decay_c*pow(delta,2));

  //we now have v, n, l, r and the material info for the sphere we hit, as well as the information about the light, so we should be able to compute
  //The calculation involves three seperate and distinct values: the ambient part, the diffuse part, and the specular part. 

  //in the given equation, the ambient part is dictated by the IaKa term, which is just from generally present light through the whole scene
  //The next term, which is (Kd(n . l)) represents the diffuse part, which is multiplied by the intensity of the light over the attenuation factor
  //the third term, Ks*(r . v)^N, is the specular part, which is also multiplied by the intensity of the light over the attenuation factor

  RGB_float color = {0, 0, 0};

  //add ambience 
  color.r += global_ambient[0]*sph->mat_ambient[0];
  color.g += global_ambient[1]*sph->mat_ambient[1];
  color.b += global_ambient[2]*sph->mat_ambient[2];

  //if the object is in the shadow of another object, stop with just ambience
  if (shadow_on && detectShadow(q, toLight, scene)) {
    return color;
  } else {

    //add diffuse
    float ndotl = vec_dot(surf_norm, toLight);
    color.r += attenuation*light1_intensity[0]*sph->mat_diffuse[0]*ndotl;
    color.g += attenuation*light1_intensity[1]*sph->mat_diffuse[1]*ndotl;
    color.b += attenuation*light1_intensity[2]*sph->mat_diffuse[2]*ndotl;

    //add specular
    color.r += attenuation*light1_intensity[0]*sph->mat_specular[0]*pow(vec_dot(reflected, toEye),sph->mat_shineness);
    color.g += attenuation*light1_intensity[1]*sph->mat_specular[1]*pow(vec_dot(reflected, toEye),sph->mat_shineness);
    color.b += attenuation*light1_intensity[2]*sph->mat_specular[2]*pow(vec_dot(reflected, toEye),sph->mat_shineness);
  }

	return color;
}


Vector genStochasticRay(Vector reflected_vector) {
	Vector temp;
	temp.x = rand();
	temp.y = rand();
	temp.z = rand();
	normalize(&temp);

	return vec_plus(reflected_vector, vec_scale(temp,0.8));
}

/************************************************************************
 * This is the recursive ray tracer - you need to implement this!
 * You should decide what arguments to use.
 ************************************************************************/
RGB_float recursive_ray_trace(Point p, Vector ray, int step) {
//
// do your thing here
//
  RGB_float color;

  if (step > 0) {
  	color.r = global_ambient[0]*0.5;
  	color.g = global_ambient[1]*0.5;
  	color.b = global_ambient[2]*0.5;
  } else {
  	color = background_clr;
  }
  RGB_float reflectedC = {0,0,0};
  RGB_float refractedC = {0,0,0};
  Point* hit = new Point;
  Point* p_hit = new Point;
  Spheres *closest_sphere;
  raysShot++;
  closest_sphere = intersect_scene(p, ray, scene, hit);

   	if (plane_on && intersect_plane(p, ray, chessboard_point, chessboard_normal, p_hit)){
		Vector eye_vec = get_vec(*p_hit, p);
	    Vector light_vec = get_vec(*p_hit, light1);
	    normalize(&light_vec);
	    normalize(&eye_vec);
	    color = getPlaneColor(*p_hit, eye_vec);

	    if (reflection_on && step < step_max){
	    	Vector reflected_vector = getReflectedRay(eye_vec, chessboard_normal); 
	 		normalize(&reflected_vector);
	 		//now trace THIS ray and add the resultant color * this objects reflectance value 
	 		step++;
	 		reflectedC = recursive_ray_trace(*p_hit, reflected_vector, step);
	 		color = clr_add(color, clr_scale(reflectedC, 0.9));
	 	}

  	}
 	if (closest_sphere != NULL) {
	    raysHit++;
	    Vector eye_vec = get_vec(*hit, p);
	    Vector surf_norm = sphere_normal(*hit, closest_sphere);
	    Vector light_vec = get_vec(*hit, light1);
	    normalize(&light_vec);
	    normalize(&surf_norm);
	    normalize(&eye_vec);

	    color = phong(*hit, eye_vec, surf_norm, closest_sphere);

	 	if (reflection_on && step < step_max){
	    	Vector reflected_vector = getReflectedRay(eye_vec, surf_norm); 
	 		normalize(&reflected_vector);
	 		//now trace THIS ray and add the resultant color * this objects reflectance value 
	 		step++;
	 		reflectedC = recursive_ray_trace(*hit, reflected_vector, step);
	 		step--;
	 		color = clr_add(color, clr_scale(reflectedC, closest_sphere->reflectance));
	 		if (stochastic_ref_on){
	 			//generate STO_RAYS random reflection rays, average their results
	 			int rays = 0;
	 			Vector stochastic_ray;
	 			RGB_float stoColor = black;
	 			RGB_float temp;
	 			while (rays < STO_RAYS){
	 				stochastic_ray = genStochasticRay(reflected_vector);
	 				normalize(&stochastic_ray);
	 				temp = recursive_ray_trace(*hit, stochastic_ray, step_max);
	 				stoColor = clr_add(stoColor, clr_scale(temp, 1/STO_RAYS));
	 				rays++;
	 			}
    			color.r += stoColor.r*closest_sphere->mat_diffuse[0];
    			color.g += stoColor.g*closest_sphere->mat_diffuse[1];
    			color.b += stoColor.b*closest_sphere->mat_diffuse[2];
	 		}
	 	}

	 	if (refraction_on && step < step_max){
	 		Point* refracted_hit = new Point;
	 		Vector refracted_vector = getRefractedRay(hit, eye_vec, refracted_hit, closest_sphere);
	 		step++;
	 		refractedC = recursive_ray_trace(*refracted_hit, refracted_vector, step);
	 		step--;
	 		color = clr_add(color, clr_scale(refractedC, 1.0));

	 	}


 	} 

	return color;
}

/*********************************************************************
 * This function traverses all the pixels and cast rays. It calls the
 * recursive ray tracer and assign return color to frame
 *
 * You should not need to change it except for the call to the recursive
 * ray tracer. Feel free to change other parts of the function however,
 * if you must.
 *********************************************************************/
void ray_trace() {
  int i, j;
  float x_grid_size = image_width / float(win_width);
  float y_grid_size = image_height / float(win_height);
  float x_start = -0.5 * image_width;
  float y_start = -0.5 * image_height;
  RGB_float ret_color = background_clr;
  RGB_float temp_color = {0,0,0};

  Point cur_pixel_pos;
  Point helper_pixel_pos;
  Vector ray[5];

  // ray is cast through center of pixel
  cur_pixel_pos.x = x_start + 0.5 * x_grid_size;
  cur_pixel_pos.y = y_start + 0.5 * y_grid_size;
  cur_pixel_pos.z = image_plane;
  normalize(&chessboard_normal);

  for (i=0; i<win_height; i++) {
    for (j=0; j<win_width; j++) {
    		helper_pixel_pos.x = cur_pixel_pos.x + 0.25*x_grid_size; //start with top right pixel
    		helper_pixel_pos.y = cur_pixel_pos.y + 0.25*y_grid_size;
    		helper_pixel_pos.z = cur_pixel_pos.z;
    		ray[1] = get_vec(eye_pos, helper_pixel_pos);

    		helper_pixel_pos.y = cur_pixel_pos.y - 0.25*y_grid_size; //bottom right
    		ray[2] = get_vec(eye_pos, helper_pixel_pos);

    		helper_pixel_pos.x = cur_pixel_pos.x - 0.25*x_grid_size; //bottom left
    		ray[3] = get_vec(eye_pos, helper_pixel_pos);

    		helper_pixel_pos.y = cur_pixel_pos.y + 0.25*y_grid_size; //top left
    		ray[4] = get_vec(eye_pos, helper_pixel_pos);
    		

    		ray[0] = get_vec(eye_pos, cur_pixel_pos);


      if (supersamp_on){	//do five pixels, average the pixel? scale each by 0.2, add them
		for (int k = 0; k < 5; k++){
			temp_color = recursive_ray_trace(eye_pos, ray[k], 0);
      		ret_color = clr_add(ret_color, clr_scale(temp_color, 0.2));
      	}
      } else {
      	ret_color = recursive_ray_trace(eye_pos, ray[0], 0);
      }

      // Parallel rays can be cast instead using below
      //
      // ray.x = ray.y = 0;
      // ray.z = -1.0;
      // ret_color = recursive_ray_trace(cur_pixel_pos, ray, 1);

      frame[i][j][0] = GLfloat(ret_color.r);
      frame[i][j][1] = GLfloat(ret_color.g);
      frame[i][j][2] = GLfloat(ret_color.b);
      ret_color = black;

      cur_pixel_pos.x += x_grid_size;
    }

    cur_pixel_pos.y += y_grid_size;
    cur_pixel_pos.x = x_start;
  }
}
