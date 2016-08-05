#include "sphere.h"
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <iostream>


/**********************************************************************
 * This function intersects a ray with a given sphere 'sph'. You should
 * use the parametric representation of a line and do the intersection.
 * The function should return the parameter value for the intersection, 
 * which will be compared with others to determine which intersection
 * is closest. The value -1.0 is returned if there is no intersection
 *
 * If there is an intersection, the point of intersection should be
 * stored in the "hit" variable
 **********************************************************************/
float intersect_sphere(Point o, Vector u, Spheres *sph, Point *hit) {
  //use equation : |o + tu - sph.center|^2 - sph.radius^2 = 0; 
  //solve for t using quadratic formula, and determine intersecting points. 
  float a, b, c, cHelper;
  Vector shiftedCenter = get_vec(sph->center, o); //(o - sph.center)
  cHelper = vec_len(shiftedCenter);

  a = vec_dot(u,u);
  b = 2*vec_dot(u, shiftedCenter);
  c = cHelper*cHelper - sph->radius*sph->radius;
  float t, t1, t2;

  //now we think about the discriminant created by these variables in the quadratic equation
  //discriminant = b^2 - 4AC; if this is positive, we have two real intersections of the sphere. If this is 0, 
  //we have 1 real intersection of the sphere. If this is negative, we have no real intersections of the sphere.

  float discriminant = b*b - 4*a*c;


  //the ray is oriented such that it would pierce the sphere. There are TWO points, and we need to figure 
  //which one is the point to return. This is easy, as it's the smaller T.
  if (discriminant>0){

    t1 = (-b + sqrt(discriminant))/(2.0*a);
    t2 = (-b - sqrt(discriminant))/(2.0*a);

    t = std::min(t1, t2); //does not handle if we are INSIDE a sphere to start
    if (t < 0.001){ 
                    //if it is negative, that means the sphere is BEHIND the vector;
      return -1.0;
    } 
    

    hit->x = o.x + t*u.x;
    hit->y = o.y + t*u.y;
    hit->z = o.z + t*u.z;
    
    return t;
    

  } else if (discriminant == 0) {  //there is one intersection and it is because the ray forms a perfect tangent of the sphere
  
    t = -b/(2*a);
    if (t < 0.001){                   //if it is negative, that means the sphere is BEHIND the vector;
      return -1.0;
    } 
    hit->x = o.x + t*u.x;
    hit->y = o.y + t*u.y;
    hit->z = o.z + t*u.z;
    
    return t;

  } else {                        //no intersection, return -1.0
    return -1.0;
  }
}

/*********************************************************************
 * This function returns a pointer to the sphere object that the
 * ray intersects first; NULL if no intersection. You should decide
 * which arguments to use for the function. For exmaple, note that you
 * should return the point of intersection to the calling function.
 **********************************************************************/
Spheres *intersect_scene(Point o, Vector u, Spheres* scene, Point *hit) {
  Point *hitTrack = new Point;
  Spheres *tracker = scene;
  Spheres *closest = NULL;
 
  float t = -1.0;
  float tTrack = 0.0;

  while (tracker != NULL){
    
    tTrack = intersect_sphere(o, u, tracker, hitTrack);
  
    if (tTrack > 0.0){ //if there was a hit THIS time
      if (t > 0.0){ //and there had been a hit before
        if (t > tTrack) { //tTrack is SMALLER, meaning closer
          closest = tracker;
          t = tTrack;
          hit->x = hitTrack->x;
          hit->y = hitTrack->y;
          hit->z = hitTrack->z;
        } //if tTrack is larger, we don't care, it's behind. 
      } else { //if there was no PREVIOUS hit, track this one
        t = tTrack;
        closest = tracker;
        hit->x = hitTrack->x;
        hit->y = hitTrack->y;
        hit->z = hitTrack->z;
      }
    }
    
    tracker = tracker->next; 
  }

	return closest;
}

/*****************************************************
 * This function is nearly identical to sphere intersection
 * we plan to return less overall information. Simply "true"
 * if the LIGHT vector at a point during phong modelling will
 * hit a sphere on the way to the light
 *****************************************************/

bool detectShadow(Point o, Vector u, Spheres* sph)
{
  Spheres* sphereTrack = sph;
  float a, b, c, cHelper, t1, t2, discriminant;
  Vector shiftedCenter; 

  while (sphereTrack){

    shiftedCenter = get_vec(sphereTrack->center, o); //(o - sphereTrack.center)
    cHelper = vec_len(shiftedCenter); //rather than use the vec_len function twice, 
                                      //to save additional calculations
    a = vec_dot(u,u);
    b = 2*vec_dot(u, shiftedCenter);
    c = cHelper*cHelper - sphereTrack->radius*sphereTrack->radius;

    discriminant = b*b - 4*a*c;
    if (discriminant > 0) { //there is one or two intersection points
      //check that t1 and t0 are not on the CURRENT circle, ie, t1 or t2 = 0.
      //Also check that t1 or t2 is positive
      t1 = (-b + sqrt(discriminant))/(2.0*a);
      t2 = (-b - sqrt(discriminant))/(2.0*a);
      if ((t1 > 0.0001) || (t2 > 0.0001)){
        return true;
      } 
    }
    sphereTrack = sphereTrack->next;
  }
  return false;
}
/*****************************************************
 * This function adds a sphere into the sphere list
 *
 * You need not change this.
 *****************************************************/
Spheres *add_sphere(Spheres *slist, Point ctr, float rad, float amb[],
		    float dif[], float spe[], float shine, 
		    float refl, float refr, int sindex) {
  Spheres *new_sphere;

  new_sphere = (Spheres *)malloc(sizeof(Spheres));
  new_sphere->index = sindex;
  new_sphere->center = ctr;
  new_sphere->radius = rad;
  (new_sphere->mat_ambient)[0] = amb[0];
  (new_sphere->mat_ambient)[1] = amb[1];
  (new_sphere->mat_ambient)[2] = amb[2];
  (new_sphere->mat_diffuse)[0] = dif[0];
  (new_sphere->mat_diffuse)[1] = dif[1];
  (new_sphere->mat_diffuse)[2] = dif[2];
  (new_sphere->mat_specular)[0] = spe[0];
  (new_sphere->mat_specular)[1] = spe[1];
  (new_sphere->mat_specular)[2] = spe[2];
  new_sphere->mat_shineness = shine;
  new_sphere->reflectance = refl;
  new_sphere->mat_refract = refr;
  new_sphere->next = NULL;

  if (slist == NULL) { // first object
    slist = new_sphere;
  } else { // insert at the beginning
    new_sphere->next = slist;
    slist = new_sphere;
  }

  return slist;
}

/******************************************
 * computes a sphere normal - done for you
 ******************************************/
Vector sphere_normal(Point q, Spheres *sph) {
  Vector rc;

  rc = get_vec(sph->center, q);
  normalize(&rc);
  return rc;
}
