#include "RayTrace.hpp"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define GLM_ENABLE_EXPERIMENTAL
#include <iomanip>
#include <algorithm> 
#include <glm/ext.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtc/epsilon.hpp>
int clampColor(int x);
string trimWhiteSpace(string str);
float alphaTri(vec3 a, vec3 b, vec3 c, vec3 direction);
vec3 calculateTransmission(float n1, float n2, vec3 direction, vec3 normal);
vec3 generateCosineWeightedPoint(float u, float v);



vec3 RayTrace::traceAndPrint(vec3 rayOrigin, vec3 ray, int depth, vec3 final_color)
{
    vec3 local_color = vec3(0,0,0);
    vec3 reflection_color = vec3(0,0,0);
    vec3 transmission_color = vec3(0,0,0);
    float fres_ref = 0;
    struct intersection inter;
    getIntersection(rayOrigin , ray, &inter, objects);

    local_color = blinnPhong(rayOrigin, ray, &inter, 1);
    vec3 absorbance = vec3(1,1,1);
    if(depth == -1)
    {
        
        cout << "  Iteration type: Primary" << endl;
        cout << "             Ray: {" << inter.rayOrigin.x << " " << inter.rayOrigin.y << " " << inter.rayOrigin.z << "} -> {" << inter.ray.x << " " << inter.ray.y << " " << inter.ray.z << "}" << endl;
        cout << " Transformed Ray: {" << inter.transformedOrigin.x << " " << inter.transformedOrigin.y << " " << inter.transformedOrigin.z << "} -> {" << inter.transformedRay.x << " " << inter.transformedRay.y << " " << inter.transformedRay.z << "}" << endl;
        cout << "      Hit Object: (ID #" << inter.objectIndex + 1 << " - " << inter.objectType << ")" << endl;
        cout << "    Intersection: {" << inter.intersectionPoint.x << " " << inter.intersectionPoint.y << " " << inter.intersectionPoint.z << "} at T = " << inter.T << endl;;
        cout << "          Normal: {" << inter.normal.x << " " << inter.normal.y << " " << inter.normal.z << "}" << endl;
        cout << "     Final Color: {" << final_color.x << " " << final_color.y << " " << final_color.z << "}" << endl;
        cout << "         Ambient: {" << inter.ambient << " " << inter.ambient << " " << inter.ambient << "}" << endl;
        cout << "         Diffuse: {" << inter.diffuseVector.x << " " << inter.diffuseVector.y << " " << inter.diffuseVector.z << "}" << endl;
        cout << "        Specular: {" << inter.specularVector.x << " " << inter.specularVector.y << " " << inter.specularVector.z << "}" << endl;
        cout << "      Reflection: {" << endl;
        cout << "      Refraction: {" ;
        cout << "   Contributions: ";
        reflection_color = traceRay(inter.intersectionPoint + inter.reflectionRay * 0.001f,inter.reflectionRay, 0, 0);
        vec3 transmission;

        //fres_ref = fresnel_reflectance(inter.ior, inter.normal, -ray);
        transmission = calculateTransmission(1, inter.ior, inter.ray, inter.normal);
        
        transmission_color = traceRay(inter.intersectionPoint + 0.001f * transmission, transmission, 0, 0);
    }
    else if(depth < 6)
    {
        if(inter.T != -1)
        {
            if(inter.filter != 0)
            {
                vec3 transmission;
                struct intersection tempinter;
                getIntersection(inter.intersectionPoint + 0.0001f*-inter.ray, -inter.ray, &tempinter, objects);
                if(tempinter.objectIndex != inter.objectIndex)
                {
                    
                    transmission = calculateTransmission(1, inter.ior, inter.ray, inter.normal);
                    
                    

                }
                else
                {
                    
                    transmission = calculateTransmission(inter.ior, 1, inter.ray, -inter.normal);
                    float distance = abs(glm::distance(tempinter.intersectionPoint, inter.intersectionPoint));

                    absorbance = (1.f - inter.pigment) * 0.15f * -distance;
                    float e = 2.71828182845904523536;
                    absorbance.x = pow(e, absorbance.x);
                    absorbance.y = pow(e, absorbance.y);
                    absorbance.z = pow(e, absorbance.z);
                    
                }
                transmission_color += traceRay(inter.intersectionPoint + 0.0001f * transmission, transmission, depth+1, 0) * absorbance;
                    
            }
            if(inter.reflection != 0)
            {
                vec3 reflectionRay;
                struct intersection tempinter;
                getIntersection(inter.intersectionPoint + 0.001f*inter.ray, inter.ray, &tempinter, objects);

                if(tempinter.objectIndex == inter.objectIndex)
                {
                    reflectionRay = ray - 2 * dot(inter.ray, -inter.normal) * -inter.normal;
                }
                else
                {
                    reflectionRay = ray - 2 * dot(inter.ray, inter.normal) * inter.normal;
                }
                
                reflection_color += traceRay(inter.intersectionPoint + 0.0001f * reflectionRay, reflectionRay, depth + 1, 0) ;
            }
        }   
        
    }
    
    
    
    float local_contribution = (1 - inter.reflection) * (1-inter.filter);
    float reflection_contribution = inter.reflection * (1-inter.filter) + inter.filter * fres_ref;
    float transmission_contribution = inter.filter * (1 - fres_ref);

    return local_color * local_contribution + reflection_color * reflection_contribution + transmission_color * transmission_contribution;
}
void RayTrace::printRays(int width, int height, int x, int y)
{
    cout << "Pixel: [" << x << ", " << y << "] Color: ";
    
    vec3 imageCoord;
    imageCoord.x = -0.5 + ((x + 0.5)/width);
    imageCoord.y = -0.5 + ((y + 0.5)/height);
    imageCoord.z = -1;
    camera.up = normalize(camera.up);
    vec3 look = camera.look_at-camera.location;
    look = normalize(look);
    vec3 direction = vec3(imageCoord.x * camera.right.x, imageCoord.y * camera.up.y, -1); 
    direction = normalize(direction);
    vec3 color = traceRay(camera.location, direction, -1, 0);
    struct intersection inter;
    getIntersection(camera.location, direction, &inter, objects);
    cout << setprecision(4);
    cout << "(" << round(color.x * 255.f) << ", " << round(color.y * 255.f) << ", " << round(color.z* 255.f) << ")" << endl;
    cout << "----" << endl;
    traceAndPrint(camera.location, direction, -1, color);
   
}
void RayTrace::pixelTrace(int width, int height, int x, int y)
{
	cout << "Pixel: [" << x << ", " << y << "] Color: ";
	
	vec3 imageCoord;
    imageCoord.x = -0.5 + ((x + 0.5)/width);
    imageCoord.y = -0.5 + ((y + 0.5)/height);
    imageCoord.z = -1;
    camera.up = normalize(camera.up);
    vec3 look = camera.look_at-camera.location;
    look = normalize(look);
   	vec3 direction = vec3(imageCoord.x * camera.right.x, imageCoord.y * camera.up.y, -1); 
   	direction = normalize(direction);


	//vec3 color = blinnPhong(camera.location, direction, 0);
vec3 color;
	cout << setprecision(4);
	cout << "(" << round(color.x * 255.f) << ", " << round(color.y * 255.f) << ", " << round(color.z* 255.f) << ")" << endl;
	cout << "o - Iteration type: Primary" << endl;
	cout << "|  Ray: {" << intersections[0].rayOrigin.x << " " << intersections[0].rayOrigin.y << " " << intersections[0].rayOrigin.z << "} -> {" << intersections[0].ray.x << " " << intersections[0].ray.y << " " << intersections[0].ray.z << "}" << endl;
	cout << "|  Hit Object ID (" << intersections[0].objectIndex + 1 << " - " << intersections[0].objectType << ") at T = ";
	cout << intersections[0].T << endl << "|  Intersection = {" << intersections[0].intersectionPoint.x << " " << intersections[0].intersectionPoint.y << " " << intersections[0].intersectionPoint.z << "}" << endl;
	cout << "|  Normal {" << intersections[0].normal.x << " " << intersections[0].normal.y << " " << intersections[0].normal.z << "}" << endl;
	cout << "|  Transformed Ray: {" << intersections[0].rayOrigin.x << " " << intersections[0].rayOrigin.y << " " << intersections[0].rayOrigin.z << "} -> {" << intersections[0].ray.x << " " << intersections[0].ray.y << " " << intersections[0].ray.z << "}" << endl;
	cout << "|  ShadowRay [0] {" << intersections[0].intersectionPoint.x << " " << intersections[0].intersectionPoint.y << " " <<  intersections[0].intersectionPoint.z << "} -> {"  << intersections[0].shadowRay.x << " " << intersections[0].shadowRay.y << " " << intersections[0].shadowRay.z << "}" << endl;
	cout << "|  Ambient: " << intersections[0].ambientVector.x << " " << intersections[0].ambientVector.y << " " << intersections[0].ambientVector.z << endl;
}

vec3 calculateTransmission(float n1, float n2, vec3 direction, vec3 normal)
{
	vec3 t;
	t = (n1 / n2) * (direction - dot(direction, normal) * normal);	
	float t3 = 1-(pow(n1/n2, 2) * (1-pow(dot(direction, normal),2)));
	t3 = sqrt(t3);
	t = t - normal * t3;
	return t;
}
void RayTrace::getIntersection(vec3 rayOrigin, vec3 direction, struct intersection *inter, vector<object> objects)
{
	if(rayOrigin == camera.location)
	{
		inter->type = "Primary";
	}
	
	vec3 closestPigment;
        	
   
    for(int i = 0; i < (signed) objects.size(); i++)
        {
        	inter->rayOrigin = rayOrigin;
			inter->ray = direction;

        	if(objects.at(i).type.compare("Sphere") == 0)
        	{
        		inter->rayOrigin = vec3(objects.at(i).invModelMatrix * vec4(rayOrigin, 1.f)) ;
        		inter->ray = vec3(objects.at(i).invModelMatrix * vec4(direction, 0.f) );
        		vec3 location = objects.at(i).location;
        		double radius = objects.at(i).length;
                int objectIndex = objects.at(i).objectIndex;
        		float a = dot(inter->ray, inter->ray);
        		float b = dot(2.f * inter->ray , inter->rayOrigin - location);
        		float c = dot(inter->rayOrigin - location, inter->rayOrigin - location) - (radius * radius);
        		float tSphere = quadraticSolver(a,b,c);
        			if(inter->T < 0 && tSphere > 0)
        			{
                        inter->transformedOrigin = inter->rayOrigin;
                        inter->transformedRay = inter->ray;
        				inter->T = tSphere;
        				inter->intersectionPoint = rayOrigin + (tSphere * direction);
        				inter->normal = (inter->rayOrigin + inter->ray * tSphere) - location;
						inter->normal = vec3(vec4(inter->normal, 0.f) * objects.at(i).invModelMatrix);
                        inter->normal = normalize(inter->normal);
        				inter->objectType = "Sphere";
        				inter->objectIndex = objectIndex;
        				vec3 reflectionRay = inter->ray - 2 * dot(inter->ray, inter->normal) * inter->normal;
        				inter->reflectionRay = reflectionRay;
                        inter->obj = &objects.at(i);
        			}
        			else
        			{
        				if(tSphere < inter->T && tSphere > 0 )
        				{
                            inter->transformedOrigin = inter->rayOrigin;
                            inter->transformedRay = inter->ray;
        					inter->T = tSphere;
        					inter->intersectionPoint = rayOrigin + (tSphere *  direction);
        					inter->normal = (inter->rayOrigin + inter->ray * tSphere) - location;
							inter->normal = vec3(vec4(inter->normal, 0.f) * objects.at(i).invModelMatrix);
                            inter->normal = normalize(inter->normal);
        					inter->objectType = "Sphere";
        					inter->objectIndex = objectIndex;
        					vec3 reflectionRay = inter->ray - 2 * dot(inter->ray, inter->normal) * inter->normal;
        					inter->reflectionRay = reflectionRay;
        					inter->obj = &objects.at(i);
        				}
        			}
        		}
        		else if(objects.at(i).type.compare("Plane") == 0)
        		{
        			inter->rayOrigin = vec3(objects.at(i).invModelMatrix * vec4(inter->rayOrigin, 1.f)); 
        			inter->ray = vec3(objects.at(i).invModelMatrix * vec4(inter->ray, 0.f));
        			vec3 location = objects.at(i).location;
        			double length = objects.at(i).length;
                    int objectIndex = objects.at(i).objectIndex;
        			float tPlane = (length - dot(inter->rayOrigin, location))/dot(inter->ray, location);
        			if(inter->T < 0 && tPlane > 0)
        			{
                        inter->transformedOrigin = inter->rayOrigin;
                        inter->transformedRay = inter->ray;
        				inter->T = tPlane;
        				inter->objectType = "Plane";
        				inter->normal = location;
        				inter->normal = vec3(vec4(inter->normal, 0.f) * objects.at(i).invModelMatrix);
                        inter->normal = normalize(inter->normal);
        				inter->objectIndex = objectIndex;
        				inter->intersectionPoint = rayOrigin + inter->T * direction- 0.001f * direction;
        				vec3 reflectionRay = inter->ray - 2 * dot(inter->ray, inter->normal) * inter->normal;
        				inter->reflectionRay = reflectionRay;
        				inter->obj = &objects.at(i);

        			}
        			else
        			{
        					if(tPlane < inter->T && tPlane > 0)
        					{
                                inter->transformedOrigin = inter->rayOrigin;
                                inter->transformedRay = inter->ray;
        						inter->T = tPlane;
        						inter->objectType = "Plane";
        						inter->normal = location;
        						inter->normal = vec3(vec4(inter->normal, 0.f) * objects.at(i).invModelMatrix);
                                inter->normal = normalize(inter->normal);
        						inter->objectIndex = objectIndex;
        						inter->intersectionPoint = rayOrigin + inter->T * direction - 0.001f * direction;
        						vec3 reflectionRay = inter->ray - 2 * dot(inter->ray, inter->normal) * inter->normal;
        						inter->reflectionRay = reflectionRay;
        						inter->obj = &objects.at(i);
        					}
        			}
        		}
        		else if(objects.at(i).type.compare("Triangle") == 0)
        		{
        			vec3 location = objects.at(i).location;
        			vec3 location2 = objects.at(i).location1;
        			vec3 location3 = objects.at(i).location2;
                    int objectIndex = objects.at(i).objectIndex;
                    vec3 A = location2 - location;
                    vec3 B = location3 - location;
                    vec3 normal = glm::cross(A, B);
                    normal = normalize(normal);
        			float aTri = alphaTri(location, location2, location3, direction);
        			float bTriangle = betaTriangle(location, location2, location3, direction, rayOrigin);
        			float gTriangle = gammaTriangle(location, location2, location3, direction, rayOrigin); 
        			float tTri = tTriangle(location, location2, location3, direction, rayOrigin);
        			bTriangle = bTriangle / aTri;
        			gTriangle = gTriangle / aTri;
                    tTri = tTri / aTri;

        			bool isHit = false;
        			if(tTri <= 0 || (tTri > inter->T && inter->T > 0))
        			{
        				isHit = true;
        			}
        			else if(gTriangle < 0 || gTriangle > 1)
        			{
        				isHit = true;
        			}
        			else if(bTriangle < 0 || bTriangle > 1 - gTriangle)
        			{
        				isHit = true;
        			}
        			if(!isHit)
        			{
        				if(inter->T < 0)
        				{

        					inter->T = tTri;
        					inter->objectType = "Triangle";
                            inter->transformedRay = inter->ray;
                            inter->transformedOrigin = inter->rayOrigin;
        					inter->normal = normal;
        					inter->objectIndex = objectIndex;
        					inter->intersectionPoint = rayOrigin + ((tTri-0.001f) * direction);
        					vec3 reflectionRay = inter->ray - 2 * dot(inter->ray, inter->normal) * inter->normal;
        					inter->reflectionRay = reflectionRay;    
                            inter->obj = &objects.at(i);    					
        				}
        				else
        				{
        					if(tTri < inter->T )
        					{
        						inter->T = tTri;
        						inter->objectType = "Triangle";
                                inter->transformedRay = inter->ray;
                                inter->transformedOrigin = inter->rayOrigin;
        						inter->normal = normal;
        						inter->objectIndex = objectIndex;
        						inter->intersectionPoint = inter->rayOrigin + ((tTri - 0.001f) * direction);
        						vec3 reflectionRay = inter->ray - 2 * dot(inter->ray, inter->normal) * inter->normal;
        						inter->reflectionRay = reflectionRay;
                                inter->obj = &objects.at(i);
        					}
        				}
        			}
        		}
                else if(objects.at(i).type.compare("Box") == 0)
                {
                    inter->rayOrigin = vec3(objects.at(i).invModelMatrix * vec4(rayOrigin, 1.f)); 
                    inter->ray = vec3(objects.at(i).invModelMatrix * vec4(direction, 0.f));
                    vec3 minimum = objects.at(i).location;
                    vec3 maximum = objects.at(i).location1;
                    int objectIndex = objects.at(i).objectIndex;
                   
                    float tBox = getTBox(inter->rayOrigin, inter->ray, minimum, maximum);
                    
                    
             
                    if(inter->T < 0 && tBox > 0)
                    {
                        inter->transformedOrigin = inter->rayOrigin;
                        inter->transformedRay = inter->ray;
                        inter->intersectionPoint = inter->rayOrigin + tBox * inter->ray;

                        inter->T = tBox;
                        inter->objectType = "Box";
                       
                        if(glm::epsilonEqual(inter->intersectionPoint.x, minimum.x, 0.001f))
                        {
                            inter->normal = vec3(-1, 0, 0);
                        }
                        else if(glm::epsilonEqual(inter->intersectionPoint.x ,maximum.x, 0.001f))
                        {
                            inter->normal = vec3(1, 0, 0);
                        }
                        else if(glm::epsilonEqual(inter->intersectionPoint.y,minimum.y, 0.001f))
                        {
                            inter->normal = vec3(0, -1, 0);
                        }
                        else if(glm::epsilonEqual(inter->intersectionPoint.y , maximum.y, 0.001f))
                        {
                            inter->normal = vec3(0, 1, 0);
                        }
                        else if(glm::epsilonEqual(inter->intersectionPoint.z,minimum.z, 0.001f))
                        {
                            inter->normal = vec3(0, 0, -1);
                        }
                        else if(glm::epsilonEqual(inter->intersectionPoint.z,maximum.z, 0.001f))
                        {
                            inter->normal = vec3(0, 0, 1);
                        }
                        inter->normal = vec3(vec4(inter->normal, 0.f) * objects.at(i).invModelMatrix);

                        inter->normal = normalize(inter->normal);
                        inter->objectIndex = objectIndex;
                        
                        vec3 reflectionRay = inter->ray - 2 * dot(inter->ray, inter->normal) * inter->normal;
                        reflectionRay = vec3(objects.at(i).invModelMatrix * vec4(reflectionRay, 0.f));
                        inter->reflectionRay = reflectionRay;
                        inter->obj = &objects.at(i);
                        inter->intersectionPoint = rayOrigin + tBox * direction;
                        

                    }
                    else
                    {
                        if(tBox < inter->T && tBox > 0)
                        {
                            inter->transformedOrigin = inter->rayOrigin;
                            inter->transformedRay = inter->ray;
                            inter->intersectionPoint = inter->rayOrigin + tBox * inter->ray;
                            inter->T = tBox;
                            inter->objectType = "Box";

                            if(glm::epsilonEqual(inter->intersectionPoint.x, minimum.x, 0.001f))
                            {
                                inter->normal = vec3(-1, 0, 0);
                            }
                            else if(glm::epsilonEqual(inter->intersectionPoint.x ,maximum.x, 0.001f))
                            {
                                inter->normal = vec3(1, 0, 0);
                            }
                            else if(glm::epsilonEqual(inter->intersectionPoint.y,minimum.y, 0.001f))
                            {
                                inter->normal = vec3(0, -1, 0);
                            }
                            else if(glm::epsilonEqual(inter->intersectionPoint.y , maximum.y, 0.001f))
                            {
                                inter->normal = vec3(0, 1, 0);
                            }
                            else if(glm::epsilonEqual(inter->intersectionPoint.z,minimum.z, 0.001f))
                            {
                                inter->normal = vec3(0, 0, -1);
                            }
                            else if(glm::epsilonEqual(inter->intersectionPoint.z,maximum.z, 0.001f))
                            {
                                inter->normal = vec3(0, 0, 1);
                            }
                            inter->normal = vec3(vec4(inter->normal, 0.f) * objects.at(i).invModelMatrix);

                         
                            inter->normal = normalize(inter->normal);
                            inter->objectIndex = i;
                            
                            inter->intersectionPoint = rayOrigin + tBox * direction;
                            
                           
                        vec3 reflectionRay = inter->ray - 2 * dot(inter->ray, inter->normal) * inter->normal;
                        reflectionRay = vec3(objects.at(i).invModelMatrix * vec4(reflectionRay, 0.f));
                        inter->reflectionRay = reflectionRay;
                        inter->obj = &objects.at(i);
                                
                                

                        }
                    }
                }
        	}
}



vec3 RayTrace::blinnPhong(vec3 rayOrigin, vec3 direction, struct intersection* inter, int type)
{
	   vec3 point = inter->intersectionPoint - 0.01f * inter->transformedRay;
	   vec3 color = vec3(0, 0, 0);

	   
       
		if(inter->T != -1)
		{
            if(type == 4 || type == 5)
            {
                int carloSize;
                if(type == 4)
                {
                    carloSize = 128;
                }
                else
                {
                    carloSize = 32;
                }
                
                vec3 carloFactor = vec3(0, 0, 0);
                vector<vec3> samplePoints = generate_hemisphere_sample_points(carloSize);
                for(int i = 0; i < (signed) samplePoints.size(); i++)
                {
                    vec3 carloRay = alignSampleVector(samplePoints.at(i), inter->normal);
                    carloFactor += traceRay(point + carloRay * 0.001f, carloRay, -1, 5) * dot(carloRay, inter->normal);
                }
                
                color = inter->obj->pigment * carloFactor * .035;
            }
            else
            {
                color = inter->obj->pigment * inter->obj->ambient;
            }

			
			
			for(int i = 0; i < (signed) light_source.size(); i++)
			{
               
				bool in_shadow = false;
				vec3 lightVec = -point + light_source.at(i).location;
				lightVec = normalize(lightVec);
				inter->shadowRay = lightVec;
				struct intersection shadow;
                if(type == 0 || type == 10)
                {
                    bvh_traversal(light_source.at(i).location, -lightVec, root, &shadow);
                }
                else
                {
                    getIntersection(light_source.at(i).location , -lightVec, &shadow, objects);
                }
				
				if(shadow.T > 0)
				{
					if(shadow.T < abs(distance(light_source.at(i).location, point)))
					{
						in_shadow = true;
						inter->inshadow = true;
                        if(type == 9)
                        {
                            if(inter->obj->filter == 0)
                            {
                                color = color + color *  dot(-shadow.normal, shadow.transformedRay) * shadow.obj->filter;
                            }
                            else
                            {
                                in_shadow = false;
                                inter->inshadow = false;
                            }
                            
                        }
                        if(type == 10)
                        {

                        }
                        
					}
				}
				if(!in_shadow)
				{
					float halfDotNormal = 0;
					vec3 Kd = inter->obj->pigment * (float) inter->obj->diffuse;
                    inter->diffuseVector = Kd;
					float normalDotLightVector = dot(lightVec, inter->normal);
                   
        			vec3 Ls = light_source.at(i).pigment;
        			vec3 Ld = light_source.at(i).pigment;
        			vec3 Ks = inter->obj->pigment * (float) inter->obj->specular;
                    
        			if(normalDotLightVector < 0)
        			{
        				normalDotLightVector = 0;
        			}
        			if(normalDotLightVector > 0)
        			{        					
        				vec3 view = inter->transformedOrigin - inter->intersectionPoint;
        				view = normalize(view);
        				vec3 half = view + lightVec;
        				half = normalize(half);
       					halfDotNormal = dot(half, inter->normal);
       					
       					float roughness = 2/pow(inter->obj->roughness,2) - 2;
       					
       					if(halfDotNormal < 0)
       					{
       						halfDotNormal = 0;
       					}
       					halfDotNormal = pow(halfDotNormal, roughness);
       				}	
        			color += ((Kd * normalDotLightVector * Ld) + (Ks * halfDotNormal * Ls));
				}
			}
		}
        return color;
}

void RayTrace::makeBVHTree(struct bvh_node* root, vector<RayTrace::object> objs, int axis)
{

    if(objs.size() <= 1)
    {
        root->obj = objs;
        calculateBBox(root, objs);
        return;
    }
    objs = sort_objects_on_axis(objs, axis + 6);
    root->left = new bvh_node;
    root->right = new bvh_node;
    makeBVHTree(root->left, left_half_arr(objs), (axis + 1) % 3 );
    makeBVHTree(root->right, right_half_arr(objs), (axis + 1) % 3);
    calculateBBox(root, objs);
}
void RayTrace::initializeBVH_Node(struct bvh_node* root)
{
    root->left = new bvh_node;
    root->right = new bvh_node;
}
void RayTrace::bvh_traversal(vec3 rayOrigin, vec3 ray, struct bvh_node* bvh, struct intersection* inter)
{
    if(getTBox(rayOrigin, ray, bvh->minimum, bvh->maximum) != -1)
    {
        if(bvh->left != NULL && getTBox(rayOrigin, ray, bvh->left->minimum, bvh->left->maximum) != -1)
        {
            bvh_traversal(rayOrigin, ray, bvh->left, inter);                  
        }
        if(bvh->right != NULL && getTBox(rayOrigin, ray, bvh->right->minimum, bvh->right->maximum) != -1)
        {
            bvh_traversal(rayOrigin, ray, bvh->right, inter);           
        }
        if(bvh->obj.size() >= 1)
        {
            getIntersection(rayOrigin, ray, inter, bvh->obj);
        }
    }
}
float RayTrace::getTBox(vec3 rayOrigin, vec3 ray, vec3 minimum, vec3 maximum)
{
    float tgmin = -std::numeric_limits<float>::infinity();
    float tgmax = std::numeric_limits<float>::infinity();

    float t1 = (minimum.x - rayOrigin.x) / ray.x;
    float t2 = (maximum.x - rayOrigin.x) / ray.x;

    if(ray.x == 0)
    {
        if(rayOrigin.x >= minimum.x || rayOrigin.x <= maximum.x)
        {
            return -1;
        }
    }
    if(t1 > t2)
    {
        swap(t1, t2);
    }

    if(t1 > tgmin)
    {
        tgmin = t1;
    }
    if(t2 < tgmax)
    {
        tgmax = t2;
    }

    float t3 = (minimum.y - rayOrigin.y) / ray.y;
    float t4 = (maximum.y - rayOrigin.y) / ray.y;

    if(ray.y == 0)
    {
        if(rayOrigin.y >= minimum.y || rayOrigin.y <= maximum.y)
        {
            return -1;
        }
    }

    if(t3 > t4)
    {
        swap(t3, t4);
    }

    if(t3 > tgmin)
    {
        tgmin = t3;
    }

    if(t4 < tgmax)
    {
        tgmax = t4;
    }

    float t5 = (minimum.z - rayOrigin.z) / ray.z;
    float t6 = (maximum.z - rayOrigin.z) / ray.z;
    if(ray.z == 0)
    {
        if(rayOrigin.z >= minimum.z || rayOrigin.z <= maximum.z)
        {
            return -1;
        }
    }

    if(t5 > t6)
    {
        swap(t5, t6);
    }

    if(t5 > tgmin)
    {
        tgmin = t5;
    }

    if(t6 < tgmax)
    {
        tgmax = t6;
    }

    if(tgmin >= tgmax)
    {
        return -1;
    }
    if(tgmax < 0)
    {
        return -1;
    }

    if(tgmin > 0 && tgmin < std::numeric_limits<float>::infinity())
    {
        return tgmin;
    }
    else if(tgmax < std::numeric_limits<float>::infinity())
    {
        return tgmax;
    }
    else
    {
        return -1;
    }
}

vec3 generateCosineWeightedPoint(float u, float v) {
    float radial = sqrt(u);
    float theta = 2.0 * 3.141592596 * v;

    float x = radial * cos(theta);
    float y = radial * sin(theta);

    return vec3(x, y, sqrt(1 - u));
}
vector<vec3> RayTrace::generate_hemisphere_sample_points(int numPoints)
{
    vector<vec3> ret;
    for(int i = 0; i < numPoints; i++)
    {
        const float u = rand() / (float) RAND_MAX;
        const float v = rand() / (float) RAND_MAX;
        vec3 point = generateCosineWeightedPoint(u, v);
        ret.push_back(point);
    }
    return ret;
}

vec3 RayTrace::alignSampleVector(vec3 sample, vec3 normal)
{
    float angle = acos(dot(camera.up, normal));
    vec3 axis = cross(camera.up, normal);

    mat4 Rotation = mat4(1.f);
    if(angle != 0)
    {
        Rotation = glm::rotate(glm::mat4(1.f), glm::radians(angle), axis) * Rotation;
        return vec3(Rotation * vec4(sample, 1.f));
    }
    return sample;
    
}
vec3 RayTrace::traceRay(vec3 rayOrigin, vec3 ray, int depth, int type)
{
	vec3 local_color = vec3(0,0,0);
	vec3 reflection_color = vec3(0,0,0);
	vec3 transmission_color = vec3(0,0,0);
    vec3 carloFactor = vec3(1, 1, 1);
 
	float fres_ref = 0;
	struct intersection inter;
    if(type == 0 || type == 10) // bvh-traversal
    {
        bvh_traversal(rayOrigin, ray, root, &inter);
        local_color = blinnPhong(inter.transformedOrigin, inter.transformedRay, &inter, 0);
    }
    if(type == 4) // carlo-ray bounce 1
    {
        getIntersection(rayOrigin , ray, &inter, objects); 
        local_color = blinnPhong(inter.transformedOrigin, inter.transformedRay, &inter, 5); 
    }
    else if(type == 5) // carlo-ray bounce 2
    {
        getIntersection(rayOrigin, ray, &inter, objects);
        local_color = blinnPhong(inter.transformedOrigin, inter.transformedRay, &inter, 6);
        depth = 6;
    }
    else if(type != 0)
    {
        getIntersection(rayOrigin , ray, &inter, objects);
        local_color = blinnPhong(inter.transformedOrigin, inter.transformedRay, &inter, type); 
    }

	vec3 absorbance = vec3(1,1,1);
	if(depth == -1)
	{
        if(inter.T != -1)
        {
            if(inter.obj->reflection != 0)
            {
                reflection_color = traceRay(inter.intersectionPoint + inter.reflectionRay * 0.001f, inter.reflectionRay, 0, type);// * inter.obj->pigment;
            }
        
            vec3 transmission;
            if(inter.obj->filter != 0 && type == 3) // fresnel-reflectance
            {
                fres_ref = fresnel_reflectance(inter.obj->ior, inter.normal, -inter.ray);
            }   
            if(inter.obj->filter != 0)
            {
                transmission = calculateTransmission(1, inter.obj->ior, inter.ray, inter.normal);
                transmission_color = traceRay(inter.intersectionPoint + 0.001f * transmission, transmission, 0, type);
            }
        }
        
	}
	else if(depth < 6)
	{
		if(inter.T != -1)
		{
			if(inter.obj->filter != 0)
			{
				vec3 transmission;
				struct intersection tempinter;
                if(type == 0 || type == 10)
                {
                    bvh_traversal(inter.intersectionPoint + 0.0001f* -inter.ray, -inter.transformedRay, root, &tempinter);
                }
                else
                {
                    getIntersection(inter.intersectionPoint + 0.0001f*-inter.ray, -inter.transformedRay, &tempinter, objects);
                }
				

        		if(tempinter.objectIndex != inter.objectIndex)
        		{
        			
           			transmission = calculateTransmission(1, inter.obj->ior, inter.transformedRay, inter.normal);
        			fres_ref = fresnel_reflectance(inter.obj->ior, inter.normal, -inter.transformedRay);
        		}
        		else
        		{
        			transmission = calculateTransmission(inter.obj->ior, 1, inter.transformedRay, -inter.normal);
        			float distance = abs(glm::distance(tempinter.intersectionPoint, inter.intersectionPoint));
                    if(type == 3)
                    {
                          fres_ref = fresnel_reflectance(inter.obj->ior, -inter.normal, -inter.ray);
                    }
                  
        			absorbance = (1.f - inter.obj->pigment) * 0.15f * -distance;
        			float e = 2.71828182845904523536;
        			absorbance.x = pow(e, absorbance.x);
        			absorbance.y = pow(e, absorbance.y);
        			absorbance.z = pow(e, absorbance.z);
        		}
				transmission_color += traceRay(inter.intersectionPoint + 0.0001f * transmission, transmission, depth+1, 1) * absorbance;
			}
			if(inter.obj->reflection != 0)
			{
				vec3 reflectionRay;
				struct intersection tempinter;
                if(type == 0 || type == 10)
                {
                    bvh_traversal(inter.intersectionPoint + 0.0001f* -inter.transformedRay, -inter.transformedRay, root, &tempinter);
                }
                else
                {
                    getIntersection(inter.intersectionPoint + 0.001f*-inter.ray, -inter.ray, &tempinter, objects);
                }
        		if(tempinter.objectIndex != inter.objectIndex)
        		{
					reflectionRay = ray - 2 * dot(inter.transformedRay, -inter.normal) * -inter.normal;
        		}
        		else
        		{
        			reflectionRay = ray - 2 * dot(inter.transformedRay, inter.normal) * inter.normal;
        		}
        		//reflectionRay = inter.reflectionRay;
				reflection_color += traceRay(inter.intersectionPoint + 0.0001f * reflectionRay, reflectionRay, depth + 1, type);
			}
		}		
	}
    if(inter.T != -1)
    {
        float local_contribution = (1 - inter.obj->reflection) * (1-inter.obj->filter) ;
        float reflection_contribution = inter.obj->reflection * (1-inter.obj->filter) + inter.obj->filter * fres_ref;
        float transmission_contribution = inter.obj->filter * (1 - fres_ref);

        return local_color * local_contribution * carloFactor + reflection_color * reflection_contribution + transmission_color  * transmission_contribution;
    }
	else
    {
        return vec3(0, 0, 0);
    }
}

float RayTrace::fresnel_reflectance(float ior, vec3 normal, vec3 ray)
{
	float F0 = pow(ior - 1, 2) / pow(ior + 1, 2);
	
	float fresnel;
	fresnel = F0 + ((1-F0) * pow(1-dot(normal, ray),5));
	return fresnel;
}
float RayTrace::checkForIntersection(vec3 point, vec3 lvector)
{
			vec3 closestPigment;
    		string type;
    		float t = -1; 
    		point = point + (0.001f * lvector);
        	for(int i = 0; i < (signed) objects.size(); i++)
        	{
        		if(objects.at(i).type.compare("Sphere") == 0)
        		{
        			vec3 location = vec3(objects.at(i).location.x, objects.at(i).location.y, objects.at(i).location.z);
        			vec3 pigment = objects.at(i).pigment;
        			double radius = objects.at(i).length;
        		
        			float a = dot(lvector, lvector);
        			float b = dot(2.f * lvector , point - location);
        			float c = dot(point - location, point - location) - (radius * radius);
        			float tSphere = quadraticSolver(a,b,c);
        			if(t == -1 && tSphere > 0)
        			{
        				t = tSphere;
        				closestPigment = pigment;
        				type = "Sphere";
        			}
        			else
        			{
        				if(tSphere < t && tSphere > 0 )
        				{

        					t = tSphere;
        					closestPigment = pigment;
        					type = "Sphere";
        				}
        			}
        		}
        		else if(objects.at(i).type.compare("Plane") == 0)
        		{
        			vec3 location = objects.at(i).location;
        			vec3 pigment = objects.at(i).pigment;
        			double length = objects.at(i).length;
        			float tPlane = (length - dot(point, location))/dot(lvector, location);
        			if(t == -1 && tPlane > 0)
        			{
        				t = tPlane;
        				closestPigment = pigment;
        				type = "Plane";
        			}
        			else
        			{
        					if(tPlane < t && tPlane > 0)
        					{
        						t = tPlane;
        						closestPigment = pigment;
        						type = "Plane";
        					}
        			}
        		}
        	}
        	
  
        	
        	return t;
        		
}
	

void RayTrace::printFirstHit(int width, int height, int x, int y)
{
	const string fileName = "output.png";
	const glm::ivec2 size = glm::ivec2(width, height);

	vec3 cameraLocation;

	
        	vec3 imageCoord;
        	imageCoord.x = -0.5 + ((x + 0.5)/size.x);
        	imageCoord.y = -0.5 + ((y + 0.5)/size.y);
        	imageCoord.z = -1;
        	camera.up = normalize(camera.up);
        	vec3 look = camera.look_at-camera.location;
        	look = normalize(look);

        	vec3 direction = vec3(imageCoord.x * camera.right.x, imageCoord.y * camera.up.y, -1); 
        	direction = normalize(direction);


        	vec3 closestPigment;
        	string type;
        	float t = -1; 
        	for(int i = 0; i < (signed) objects.size(); i++)
        	{
        		if(objects.at(i).type.compare("Sphere") == 0)
        		{
        			vec3 location = vec3(objects.at(i).location.x, objects.at(i).location.y, objects.at(i).location.z);
        			vec3 pigment = objects.at(i).pigment;
        			double radius = objects.at(i).length;
        		
        			float a = dot(direction, direction);
        			float b = dot(2.f * direction , camera.location - location);
        			float c = dot(camera.location - location, camera.location - location) - (radius * radius);
        			float tSphere = quadraticSolver(a,b,c);
     
        			if(t == -1 && tSphere > 0)
        			{
        				t = tSphere;
        				closestPigment = pigment;
        				type = "Sphere";
        			}
        			else
        			{
        				if(tSphere < t && tSphere > 0 )
        				{

        					t = tSphere;
        					closestPigment = pigment;
        					type = "Sphere";
        					

        				}
        			}
        		}
        		else if(objects.at(i).type.compare("Plane") == 0)
        		{
        			vec3 location = objects.at(i).location;
        			vec3 pigment = objects.at(i).pigment;
        			double length = objects.at(i).length;
        			float tPlane = (length - dot(camera.location, location))/dot(direction, location);
        			if(t == -1 && tPlane > 0)
        			{
        				t = tPlane;
        				closestPigment = pigment;
        				type = "Plane";
        			}
        			else
        			{
        					if(tPlane < t && tPlane > 0)
        					{
        						t = tPlane;
        						closestPigment = pigment;
        						type = "Plane";
        					}
        			}
        		}
        	}
        cout << setprecision(4);
    	cout << "Pixel: [" << x << ", " << y << "] Ray: {" << camera.location.x << " " << camera.location.y << " "
    	<< camera.location.z << "} -> {" << direction.x << " " << direction.y << " " << direction.z << "}" << endl;
        	if(t == -1)
        	{
        		cout << "No Hit" << endl;
        	}
        	else
        	{
        		cout << "T = " << t << endl << "Object Type: " << type << endl;
        		cout << "Color: " << closestPigment.x << " " << closestPigment.y << " " << closestPigment.z << endl;
        	}
	
}
void RayTrace::printSceneInfo()
{
	int j = 0;
	cout << setprecision(4);
	cout << "Camera:" << endl;
	cout << "- Location: {" << camera.location.x << " " << camera.location.y << " " << camera.location.z << "}" << endl;
	cout << "- Up: {" << camera.up.x << " " << camera.up.y << " " << camera.up.z << "}" << endl;
	cout << "- Right: {" << camera.right.x << " " << camera.right.y << " " << camera.right.z << "}" << endl;
	cout << "- Look at: {" << camera.look_at.x << " " << camera.look_at.y << " " << camera.look_at.z << "}" << endl;
	cout << endl << "---" << endl << endl;
	cout << light_source.size() << " light(s)" << endl << endl;
	for(int i = 0; i < (signed) light_source.size(); i ++)
	{
		cout << "Light[" << i << "]:" << endl;
		cout << "- Location: {" << light_source.at(i).location.x << " " << light_source.at(i).location.y << " " << light_source.at(i).location.z << "}" << endl;
		cout << "- Color: {" << light_source.at(i).pigment.x << " " << light_source.at(i).pigment.y << " " << light_source.at(i).pigment.z << "}" << endl;
		cout << endl;
	}
	cout << "---" << endl << endl;
	cout << objects.size() << " " << "object(s)";
	for(j = 0; j < (signed) objects.size(); j++)
	{
		cout << endl << endl;
		cout << "Object[" << j << "]:" << endl;
		cout << "- Type: " << objects.at(j).type << endl;
		if(objects.at(j).type.compare("Plane") == 0)
		{
			cout << "- Normal: {" << objects.at(j).location.x << " " << objects.at(j).location.y << " " << objects.at(j).location.z << "}" << endl;
			cout << "- Distance: " << objects.at(j).length << endl;
		}
		else if(objects.at(j).type.compare("Sphere") == 0)
		{
			cout << "- Center: {" << objects.at(j).location.x << " " << objects.at(j).location.y << " " << objects.at(j).location.z << "}" << endl;
			cout << "- Radius: " << objects.at(j).length << endl;
		}
		
		cout << "- Color: {" << objects.at(j).pigment.x << " " << objects.at(j).pigment.y << " " << objects.at(j).pigment.z << "}" << endl;
		cout << "- Material:" << endl;
		cout << "  - Ambient: " << objects.at(j).ambient << endl;
		cout << "  - Diffuse: " << objects.at(j).diffuse;
	}
	
	cout << endl;
}
void RayTrace::printRayInfo(int width, int height, int x, int y)
{
	const glm::ivec2 size = glm::ivec2(width, height);
	vec3 imageCoord;
    imageCoord.x = -0.5 + ((x + 0.5)/size.x);
    imageCoord.y = -0.5 + ((y + 0.5)/size.y);
    imageCoord.z = -1;
    camera.up = normalize(camera.up);
    vec3 look = camera.look_at - camera.location;
    look = normalize(look);
    look = -look;
    vec3 direction = vec3(imageCoord.x * camera.right) + vec3(imageCoord.y * camera.up) - look;
    direction = normalize(direction);
    cout << setprecision(4);
    cout << "Pixel: [" << x << ", " << y << "] Ray: {" << camera.location.x << " " << camera.location.y << " "
    << camera.location.z << "} -> {" << direction.x << " " << direction.y << " " << direction.z << "}" << endl;
}

vec3 RayTrace::getImageCoordinate(int x, int y, int width, int height, float subpixelM, float subpixelN, int ss)
{
	vec3 imageCoord;
	imageCoord.x = -0.5 + ((x + 0.5 + (-0.5 + ((subpixelM + 0.5)/ss)))/width);
    imageCoord.y = -0.5 + ((y + 0.5 + (-0.5 + ((subpixelN + 0.5)/ss)))/height);
    imageCoord.z = -1;
    return imageCoord;
}



void RayTrace::drawScene(int width, int height, string typ, int ss, int arg)
{
	const int numChannels = 3;
	const string fileName = "output.png";
	const glm::ivec2 size = glm::ivec2(width, height);

	vec3 cameraLocation;
    
	unsigned char *data = new unsigned char[size.x * size.y * numChannels];
	
	for (int y = 0; y < size.y; ++ y)
	{
    	for (int x = 0; x < size.x; ++ x)
   		{
   			
   			vec3 color = vec3(0,0,0);
        	camera.up = normalize(camera.up);
        	vec3 look = camera.look_at-camera.location;
        	look = normalize(look);

            for(int i = 0; i < ss ; i++)
            {

                vec3 imageCoord = getImageCoordinate(x, y, size.x, size.y,(float) i, (float) i, ss);
                vec3 direction = (imageCoord.x * camera.right) + (imageCoord.y * camera.up) + (imageCoord.z * -look);
              
                direction = normalize(direction);
                color += traceRay(camera.location, direction, -1, arg);
            }
        	color = color / (float) ss;
        	

        	int red = round(color.x*255.f);
        	int green = round(color.y*255.f);
        	int blue = round(color.z*255.f);
        	red = clampColor(red);
        	blue = clampColor(blue);
        	green = clampColor(green);
        	data[(size.x * numChannels) * (size.y - 1 - y) + numChannels * x + 0] = red;
        	data[(size.x * numChannels) * (size.y - 1 - y) + numChannels * x + 1] = green;
        	data[(size.x * numChannels) * (size.y - 1 - y) + numChannels * x + 2] = blue;

        }
    }

	stbi_write_png(fileName.c_str(), size.x, size.y, numChannels, data, size.x * numChannels);
	delete[] data;
}


        	

float RayTrace::betaTriangle(vec3 a, vec3 b, vec3 c, vec3 direction, vec3 rayOrigin)
{
	glm::mat3 t;


	t[0][0] = a.x - rayOrigin.x;
	t[0][1] = a.x - c.x;
	t[0][2] = direction.x;
	t[1][0] = a.y - rayOrigin.y;
	t[1][1] = a.y - c.y;
	t[1][2] = direction.y;
	t[2][0] = a.z - rayOrigin.z;
	t[2][1] = a.z - c.z;
	t[2][2] = direction.z;
	

	float det = determinant(t);
 
	return det;
}
float RayTrace::tTriangle(vec3 a, vec3 b, vec3 c, vec3 direction, vec3 rayOrigin)
{
	glm::mat3 t;

	t[0][0] = a.x - b.x;
	t[0][1] = a.x - c.x;
	t[0][2] = a.x - rayOrigin.x;
	t[1][0] = a.y - b.y;
	t[1][1] = a.y - c.y;
	t[1][2] = a.y - rayOrigin.y;
	t[2][0] = a.z - b.z;
	t[2][1] = a.z - c.z;
	t[2][2] = a.z - rayOrigin.z;
	
	float det = determinant(t);


  
	return det;
}

float RayTrace::gammaTriangle(vec3 a, vec3 b, vec3 c, vec3 direction, vec3 rayOrigin)
{
	glm::mat3 t;

	t[0][0] = a.x - b.x;	
	t[0][1] = a.x - rayOrigin.x;
	t[0][2] = direction.x;
	t[1][0] = a.y - b.y;
	t[1][1] = a.y - rayOrigin.y;
	t[1][2] = direction.y;
	t[2][0] = a.z - b.z;
	t[2][1] = a.z - rayOrigin.z;
	t[2][2] = direction.z;
	float det = determinant(t);
	
    
	return det;
}
float alphaTri(vec3 a, vec3 b, vec3 c, vec3 direction)
{
    glm::mat3 t;
	t[0][0] = a.x - b.x;   
    t[0][1] = a.x - c.x;
    t[0][2] = direction.x;
    t[1][0] = a.y - b.y;
    t[1][1] = a.y - c.y;
    t[1][2] = direction.y;
    t[2][0] = a.z - b.z;
    t[2][1] = a.z - c.z;
    t[2][2] = direction.z;
    float det = determinant(t);
    
	return det;
}
int clampColor(int x)
{
	if(x < 0)
	{
		x = 0;
	}
	if(x > 255)
	{
		x = 255;
	}
	return x;
}


float RayTrace::quadraticSolver(float a, float b, float c)
{
	double D = sqrt((b*b) - (4*a*c) );
	
    float tA = ( -b + D)/(2*a);
    float tB = ( -b - D)/(2*a);

    if(D < 0)
    {
    	return -1;
    }
    if(tA < tB && tA > 0)
    {
    	return tA;
    }
    else if(tB < tA && tB > 0)
    {
    	return tB;
    }
    else if(tA > tB && tB < 0 && tA > 0)
    {
    	return tA;
    }
    else if(tB > tA && tA < 0 && tB > 0)
    {
    	return tB;
    }

    return -1;
}
int RayTrace::parseFile(string filename)
{
	ifstream myfile;
	string line;
	myfile.open(filename.c_str());
	if(!myfile)
	{
		cout << "File Cannot Be Found or Opened." << endl;
		return -1;
	}

    ifstream FileHandle(filename);
	string String;

	FileHandle.seekg(0, std::ios::end);
	String.reserve((uint) FileHandle.tellg());
	FileHandle.seekg(0, std::ios::beg);

	String.assign((std::istreambuf_iterator<char>(FileHandle)), std::istreambuf_iterator<char>());
	vector <string> split;
	int i = 0;

	while(i < (signed) String.length())
	{
		
		if(String.at(i) == '/' && String.at(i+1) == '/')
		{
			i = i + 2;
			while(String.at(i) != '\n')
			{
				string temp;
				temp.push_back(String.at(i));
				split.push_back(temp);
				i++;

			}
		}

		else if(String.at(i) == '{')
		{
			string type;
			while(String.at(i) != '\n')
			{
				i--;
                if(i < 0)
                {
                    break;
                }
			}
			i++;
			while(String.at(i) != '{')
			{
				type.push_back(String.at(i));
				i++;
			}
			type = trimWhiteSpace(type);
			if(type.compare("sphere") == 0)
			{
				i = parseObject(String, i, "Sphere");
			}
			if(type.compare("camera") == 0)
			{
				i = parseCamera(String, i);
			}
			if(type.compare("light_source") == 0)
			{
				i = parseLightSource(String, i);
			}
			if(type.compare("plane") == 0)
			{
				i = parseObject(String, i, "Plane");
			}
			if(type.compare("triangle") == 0)
			{
				i = parseObject(String, i, "Triangle");
			}
            if(type.compare("box") == 0)
            {
                i = parseObject(String, i, "Box");
            }
		}
		i++;

	}
    root = new bvh_node();
    
    makeBVHTree(root, objects, 0);

    return 0;
}
glm::vec3 RayTrace::parseVector(string str, double *filter)
{
	int i = 0;
	string tempX, tempY, tempZ, tempFilter;

	float x, y, z;
	if(str.at(i) != '<' || str.at(i) != '{')
	{
		i++;
	}
	while(str.at(i) != ',')
	{
		tempX.push_back(str.at(i));
		i++;
	}
	x = stof(tempX);
	i++;
	while(str.at(i) != ',')
	{
		tempY.push_back(str.at(i));
		i++;
	}
	y = stof(tempY);
	i++;
	while(str.at(i) != '>')
	{
		if(str.at(i) == ',')
		{
			i++;
			while(str.at(i) != '>')
			{
				tempFilter.push_back(str.at(i));
				i++;
			}
			*filter = stod(tempFilter);

		}
		if(str.at(i) != '>')
		{
			tempZ.push_back(str.at(i));
			i++;
		}
	}
	z = stof(tempZ);
	i++;
	return vec3(x,y,z);

}


int RayTrace::parseLightSource(const string str, int index)
{
	vec3 location;
	string loc;
	vec3 pigment;
	string pig;
	while(str.at(index) != '<')
	{
		index++;
	}
	while(str.at(index) != '>')
	{
		loc.push_back(str.at(index));
		index++;
	}
	loc.push_back('>');
	location = parseVector(loc);
	while(str.at(index) != '<')
	{
		index++;
	}
	while(str.at(index) != '>')
	{
		pig.push_back(str.at(index));
		index++;
	}
	pig.push_back('>');
	pigment = parseVector(pig);
	createLightSource(location, pigment);
	while(str.at(index) != '}')
	{
		index++;
	}
	return index;
}

int RayTrace::parseCamera(string str, int index)
{
	vec3 location, up, right, look_at;
    vec3 direction;
	string loc, u, r, la;
	while(str.at(index) != '}')
	{
		if(str.at(index) == '<')
		{
			while(str.at(index) != '\n' && str.at(index) != '{')
			{
				index--;
			}
			index++;
			string type;
			while(str.at(index) != '<')
			{
				type.push_back(str.at(index));
				index++;
			}
			type = trimWhiteSpace(type);
			if(type.compare("location") == 0)
			{
				while(str.at(index) != '<')
				{
					index++;
				}
				string temp;
				while(str.at(index) != '>')
				{
					temp.push_back(str.at(index));
					index++;
				}
				temp.push_back('>');

				location = parseVector(temp);
			}
			else if(type.compare("up") == 0)
			{
				string temp;
				while(str.at(index) != '>')
				{
					temp.push_back(str.at(index));
					index++;
				}
				temp.push_back('>');
				
				up = parseVector(temp);
			}
			else if(type.compare("right") == 0)
			{
				while(str.at(index) != '<')
				{
					index++;
				}
				string temp;
				while(str.at(index) != '>')
				{
					temp.push_back(str.at(index));
					index++;
				}
				temp.push_back('>');

				right = parseVector(temp);

			}
			else if(type.compare("look_at") == 0)
			{
				string temp;
				while(str.at(index) != '>')
				{
					temp.push_back(str.at(index));
					index++;
				}
				temp.push_back('>');
				look_at = parseVector(temp);
			}
            else if(type.compare("direction") == 0)
            {
                string temp;
                while(str.at(index) != '>')
                {
                    temp.push_back(str.at(index));
                    index++;
                }
                temp.push_back('>');
                direction = parseVector(temp);

            }
			else //replace later with finish etc etc
			{
				while(str.at(index) != '}')
				{
					index++;
				}
				index++;
			}

		}
		index++;
	}
	createCamera(location, up, right, look_at);
	return index;
}

int RayTrace::parseObject(string str, int index, string type)
{
	vec3 location;
	vec3 location1;
	vec3 location2;
	string loc;
	vec3 pigment;
	string order;
	vec3 scale = vec3(1, 1, 1);
	vec3 rotate = vec3(0, 0, 0);
	vec3 translate = vec3 (0, 0, 0);

	string rad;
	float radius;
	double ambient = 0;
	double specular = 0;
	double roughness = 0;
	double diffuse = 0;
	double reflection = 0;
	float refraction = 0;
	float ior = 0;
	double filter = 0;
	while(str.at(index) != '<')
	{
		index++;
	}
	while(str.at(index) != '>')
	{
		loc.push_back(str.at(index));
		index++;
	}
	loc.push_back('>');
	location = parseVector(loc);
	if(type.compare("Triangle")==0)
	{
		loc.clear();
		while(str.at(index) != '<')
		{
			index++;
		}
		while(str.at(index) != '>')
		{
			loc.push_back(str.at(index));
			index++;
		}
		loc.push_back('>');
		location1 = parseVector(loc);
		loc.clear();
		while(str.at(index) != '<')
		{
			index++;
		}
		while(str.at(index) != '>')
		{
			loc.push_back(str.at(index));
			index++;
		}
		loc.push_back('>');
		location2 = parseVector(loc);
	}
	else if(type.compare("Box") == 0)
    {
        loc.clear();
        while(str.at(index) != '<')
        {
            index++;
        }
        while(str.at(index) != '>')
        {
            loc.push_back(str.at(index));
            index++;
        }
        loc.push_back('>');
        location1 = parseVector(loc);
    }
    else
	{
		while(!isdigit(str.at(index)))
		{
			char temp = str.at(index);
			if (temp == '-')
			{
				break;
			}
			index++;
		}
	
		while(str.at(index) != ' ' && str.at(index) != '\n')
		{
			rad.push_back(str.at(index));
			index++;
		}
		radius = stof(rad);
	}
	
	
	while(str.at(index) != '}')
	{
		if(str.at(index) == '{' || str.at(index) == '<')
		{
			while(str.at(index) != '\n')
			{
				index--;
			}
			index++;
			string type;
			while(str.at(index) != '{' && str.at(index) != '<')
			{
				type.push_back(str.at(index));
				index++;
			}
			type = trimWhiteSpace(type);
			if(type.compare("pigment") == 0)
			{
				while(str.at(index) != '<')
				{
					index++;
				}
				string temp;
				while(str.at(index) != '>')
				{
					temp.push_back(str.at(index));
					index++;
				}
				temp.push_back('>');
				pigment = parseVector(temp, &filter);
				while(str.at(index) != '}')
				{
					index++;
				}
				index++;
			}
			else if(type.compare("finish") == 0)
			{
				index++;
				while(str.at(index) != '}')
				{
					string t;
					while(!isdigit(str.at(index)) && str.at(index) != '}')
					{
						t.push_back(str.at(index));
						index++;
					}
					t = trimWhiteSpace(t);
					if(t.compare("ambient") == 0)
					{
						string val;
						while(isdigit(str.at(index)) || str.at(index) == '.')
						{
							val.push_back(str.at(index));
							index++;
						}
						ambient = stod(val);
					}
					if(t.compare("diffuse") == 0)
					{
						string val;
						while(isdigit(str.at(index)) || str.at(index) == '.')
						{
							val.push_back(str.at(index));
							index++;
						}
						diffuse = stod(val);

					}

					if(t.compare("specular") == 0)
					{
						string val;
						while(isdigit(str.at(index)) || str.at(index) == '.')
						{
							val.push_back(str.at(index));
							index++;
						}
						specular = stod(val);
					}
					if(t.compare("roughness") == 0)
					{
						string val;
						while(isdigit(str.at(index)) || str.at(index) == '.')
						{
							val.push_back(str.at(index));
							index++;
						}
						roughness = stod(val);
					}
					if(t.compare("reflection") == 0)
					{
						string val;
						while(isdigit(str.at(index)) || str.at(index) == '.')
						{
							val.push_back(str.at(index));
							index++;
						}
						reflection = stod(val);
					}
					if(t.compare("refraction") == 0)
					{
						string val;
						while(isdigit(str.at(index)) || str.at(index) == '.')
						{
							val.push_back(str.at(index));
							index++;
						}
						refraction = stof(val);
					}
					if(t.compare("ior") == 0)
					{
						string val;
						while(isdigit(str.at(index)) || str.at(index) == '.')
						{
							val.push_back(str.at(index));
							index++;
						}
						ior = stof(val);
					}


				}


			}
			else if(type.compare("scale") == 0)
			{
				string temp;
				while(str.at(index) != '>')
				{
					temp.push_back(str.at(index));
					index++;
				}
				temp.push_back('>');
				scale = parseVector(temp);
				order.push_back('0');
				
				index++;
			}
			else if(type.compare("rotate") == 0)
			{
				string temp;
				while(str.at(index) != '>')
				{
					temp.push_back(str.at(index));
					index++;
				}
				temp.push_back('>');
				rotate = parseVector(temp);
				order.push_back('1');
				
				index++;
			}
			else if(type.compare("translate") == 0)
			{
				string temp;
				while(str.at(index) != '>')
				{
					temp.push_back(str.at(index));
					index++;
				}
				temp.push_back('>');
				order.push_back('2');
				translate = parseVector(temp);
				
				//index++;
			}
			else
			{
				while(str.at(index) != '}')
				{
					index++;
				}
				index++;
			}

		}
		index++;
	}
	if(type.compare("Triangle") == 0)
	{
		createTriangle(location, location1, location2, pigment, radius, ambient, diffuse, specular, roughness, type,filter);
	}
    else if(type.compare("Box") == 0)
    {
        createBox(location, location1, pigment, ambient, diffuse, specular, roughness, reflection, refraction, ior, type, filter, scale, rotate, translate, order );
    }
	else
	{
		createObject(location, pigment, radius, ambient, diffuse, specular, roughness, reflection, refraction, ior, type, filter, scale, rotate, translate, order);
	}
	
	return index;
	
}


string trimWhiteSpace(string str)
{
	string ret;
	for(int i = 0; i < (signed) str.length(); i++)
	{
		if(str.at(i) != ' ')
		{
			ret.push_back(str.at(i));
		}
	}
	return ret;
}

void RayTrace::createTriangle(glm::vec3 location, glm::vec3 location1, glm::vec3 location2, glm::vec3 pigment, float radius, double ambient, double diffuse, double specular, double roughness, string type, double filter)
{
	struct object sp;
	sp.location = location;
	sp.location1 = location1;
	sp.location2 = location2;
	sp.pigment = pigment;
	sp.length = radius;
	sp.ambient = ambient;
	sp.diffuse = diffuse;
	sp.type = type;
	sp.specular = specular;
	sp.roughness = roughness;
	sp.filter = filter;
    
    sp.minimum = vec3(std::min(std::min(location.x, location1.x), location2.x), std::min(std::min(location.y, location1.y), location2.y), std::max(std::max(location.z, location1.z), location2.z));
    sp.maximum = vec3(std::max(std::max(location.x, location1.x), location2.x), std::max(std::max(location.y, location1.y), location2.y), std::min(std::min(location.z, location1.z), location2.z));
   
    
    sp.min = sp.minimum;
    sp.max = sp.maximum;
    sp.min = sp.min - vec3(0.001f, 0.001f, 0.001f);
    sp.max = sp.max + vec3(0.001f, 0.001f, 0.001f);
    sp.center = (sp.min + sp.max )/ 2;
    addObject(sp);

}



void RayTrace::createBox(glm::vec3 location, glm::vec3 location1, glm::vec3 pigment, double ambient,double diffuse, double specular, double roughness, double reflection, double refraction, double ior, string type, double filter, vec3 sc, vec3 rot, vec3 trans, string order)
{
    struct object sp;
    sp.location = location;
    sp.location1 = location1;
    sp.pigment = pigment;
    sp.ambient = ambient;
    sp.diffuse = diffuse;
    sp.type = type;
    sp.specular = specular;
    sp.roughness = roughness;
    sp.reflection = reflection;
    sp.refraction = refraction;
    sp.ior = ior;
    sp.filter = filter;
    sp.scale = sc;
    sp.rotate = rot;
    sp.translate = trans;
    sp.order = order;
    mat4 ModelMatrix = mat4(1.f);
    for (int i = 0; i < (signed) order.length(); i++)
    {
        if(order.at(i) == '0')
        {
            ModelMatrix = scale(mat4(1.f), sc) * ModelMatrix;
  
        }
        else if(order.at(i) == '1')
        {
            mat4 Rotation = mat4(1.f);
            Rotation = glm::rotate(glm::mat4(1.f), glm::radians(rot.z), glm::vec3(0, 0, 1)) * Rotation;
            Rotation = glm::rotate(glm::mat4(1.f), glm::radians(rot.y), glm::vec3(0, 1, 0)) * Rotation;
            Rotation = glm::rotate(glm::mat4(1.f), glm::radians(rot.x), glm::vec3(1, 0, 0)) * Rotation;
            ModelMatrix = Rotation * ModelMatrix;

        }
        else if(order.at(i) == '2')
        {
            ModelMatrix = translate(mat4(1.f), trans) * ModelMatrix;

        }
    }

    sp.ModelMatrix = ModelMatrix;
    sp.invModelMatrix = inverse(ModelMatrix);
    sp.minimum = location;
    sp.maximum = location1;
    calculateTransformedBoundingBox(&sp);
    sp.center = (sp.minimum + sp.maximum) / 2.f;    
    sp.center = vec3(sp.ModelMatrix * vec4(sp.center, 1.f));
    addObject(sp);
}
void RayTrace::createObject(glm::vec3 location, glm::vec3 pigment, float radius, double ambient, double diffuse,double specular, double roughness, double reflection, double refraction, double ior, string type, double filter, vec3 sc, vec3 rot, vec3 trans, string order)
{
	struct object sp;
	sp.location = location;
	sp.pigment = pigment;
	sp.length = radius;
	sp.ambient = ambient;
	sp.diffuse = diffuse;
	sp.type = type;
	sp.specular = specular;
	sp.roughness = roughness;
	sp.reflection = reflection;
	sp.refraction = refraction;
	sp.ior = ior;
	sp.filter = filter;
	sp.scale = sc;
	sp.rotate = rot;
	sp.translate = trans;
	sp.order = order;
	mat4 ModelMatrix = mat4(1.f);
	for (int i = 0; i < (signed) order.length(); i++)
	{
		if(order.at(i) == '0')
		{
			ModelMatrix = scale(mat4(1.f), sc) * ModelMatrix;
  
		}
		else if(order.at(i) == '1')
		{
			mat4 Rotation = mat4(1.f);
			Rotation = glm::rotate(glm::mat4(1.f), glm::radians(rot.z), glm::vec3(0, 0, 1)) * Rotation;
			Rotation = glm::rotate(glm::mat4(1.f), glm::radians(rot.y), glm::vec3(0, 1, 0)) * Rotation;
			Rotation = glm::rotate(glm::mat4(1.f), glm::radians(rot.x), glm::vec3(1, 0, 0)) * Rotation;
			ModelMatrix = Rotation * ModelMatrix;

		}
		else if(order.at(i) == '2')
		{
			ModelMatrix = translate(mat4(1.f), trans) * ModelMatrix;

		}
	}

	sp.ModelMatrix = ModelMatrix;
	sp.invModelMatrix = inverse(ModelMatrix);
    if(sp.type.compare("Sphere") == 0)
    {
        float length = radius;
        sp.minimum = vec3(location.x - length, location.y - length, location.z + length);
        sp.maximum = vec3(location.x + length, location.y + length, location.z - length);
        sp.center = location;
        sp.center = vec3(sp.ModelMatrix * vec4(sp.center, 1.f));
        calculateTransformedBoundingBox(&sp);
        
       

    }
    if(sp.type.compare("Plane") == 0)
    {
        sp.minimum = vec3(std::numeric_limits<float>::infinity(),std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity());
        sp.maximum = vec3(-std::numeric_limits<float>::infinity(),  -std::numeric_limits<float>::infinity(),  +std::numeric_limits<float>::infinity());
        sp.min = sp.minimum;
        sp.max = sp.maximum;
        sp.center = sp.location * (float) sp.length;
    }

	addObject(sp);
}

void RayTrace::calculateTransformedBoundingBox(struct object* sp)
{
    if(sp->type.compare("Sphere") == 0)
    {
    vector<vec3> vertices;
    vector<vec3> transformedVertices;
    vec3 minimum;
    vec3 maximum;
    vertices.push_back(vec3(sp->location.x - sp->length, sp->location.x - sp->length, sp->location.y + sp->length));
    vertices.push_back(vec3(sp->location.x - sp->length, sp->location.x - sp->length, sp->location.y - sp->length));
    vertices.push_back(vec3(sp->location.x + sp->length, sp->location.x - sp->length, sp->location.y - sp->length));
    vertices.push_back(vec3(sp->location.x + sp->length, sp->location.x - sp->length, sp->location.y + sp->length));
    vertices.push_back(vec3(sp->location.x - sp->length, sp->location.x + sp->length, sp->location.y + sp->length));
    vertices.push_back(vec3(sp->location.x - sp->length, sp->location.x + sp->length, sp->location.y - sp->length));
    vertices.push_back(vec3(sp->location.x + sp->length, sp->location.x + sp->length, sp->location.y - sp->length));
    vertices.push_back(vec3(sp->location.x + sp->length, sp->location.x + sp->length, sp->location.y + sp->length));

    for(int i = 0; i < 8; i++)
    {
        transformedVertices.push_back(vec3(sp->ModelMatrix * vec4(vertices.at(i), 1.f)));

    }

    vector<vec3> sortedByX;
    sortedByX = sortVertices(transformedVertices, 0);
    minimum.x = sortedByX.at(0).x;
    maximum.x = sortedByX.at(7).x;
    
    vector<vec3> sortedByY;
    sortedByY = sortVertices(transformedVertices, 1);
    minimum.y = sortedByY.at(0).y;
    maximum.y = sortedByY.at(7).y;

    vector<vec3> sortedByZ;
    sortedByZ = sortVertices(transformedVertices, 2);
    minimum.z = sortedByZ.at(7).z;
    maximum.z = sortedByZ.at(0).z;
    
    sp->min = minimum;
    sp->max = maximum;
    }
    else if(sp->type.compare("Box") == 0)
    {
        vector<vec3> vertices;

    vector<vec3> transformedVertices;
    vec3 minimum;
    vec3 maximum;
    vertices.push_back(sp->minimum);
    vertices.push_back(vec3(sp->minimum.x, sp->minimum.y, sp->maximum.z));
    vertices.push_back(vec3(sp->maximum.x, sp->minimum.y, sp->maximum.z));
    vertices.push_back(vec3(sp->maximum.x, sp->minimum.y, sp->minimum.z));
    vertices.push_back(vec3(sp->maximum.x, sp->maximum.y, sp->minimum.z));
    vertices.push_back(vec3(sp->minimum.x, sp->maximum.y, sp->minimum.z));
    vertices.push_back(vec3(sp->minimum.x, sp->maximum.y, sp->maximum.z));
    vertices.push_back(vec3(sp->maximum));
    for(int i = 0; i < 8; i++)
    {
        transformedVertices.push_back(vec3(sp->ModelMatrix * vec4(vertices.at(i), 1.f)));

    }

    vector<vec3> sortedByX;
    sortedByX = sortVertices(transformedVertices, 0);
    minimum.x = sortedByX.at(0).x;
    maximum.x = sortedByX.at(7).x;
    
    vector<vec3> sortedByY;
    sortedByY = sortVertices(transformedVertices, 1);
    minimum.y = sortedByY.at(0).y;
    maximum.y = sortedByY.at(7).y;

    vector<vec3> sortedByZ;
    sortedByZ = sortVertices(transformedVertices, 2);
    minimum.z = sortedByZ.at(7).z;
    maximum.z = sortedByZ.at(0).z;
    
    sp->min = minimum;
    sp->max = maximum;
    }


}





void RayTrace::createCamera(glm::vec3 location, glm::vec3 up, glm::vec3 right, glm::vec3 look_at)
{
	struct camera cam;
	cam.location = location;
	cam.up = up;
	cam.right = right;
	cam.look_at = look_at;
	editCamera(cam);
}

void RayTrace::createLightSource(glm::vec3 location, glm::vec3 pigment)
{
	struct light_source ls;
	ls.location = location;
	ls.pigment = pigment;
	addLightSource(ls);
}


/*void RayTrace::printpixColorInfo(int width,int height,int  x, int y, string typ)
{
	const string fileName = "output.png";
	const glm::ivec2 size = glm::ivec2(width, height);

	vec3 cameraLocation;

	
        	vec3 imageCoord;
        	imageCoord.x = -0.5 + ((x + 0.5)/size.x);
        	imageCoord.y = -0.5 + ((y + 0.5)/size.y);
        	imageCoord.z = -1;
        	camera.up = normalize(camera.up);
        	vec3 look = camera.look_at-camera.location;
        	look = normalize(look);

        	vec3 direction = vec3(imageCoord.x * camera.right.x, imageCoord.y * camera.up.y, -1); 
        	direction = normalize(direction);


        	vec3 closestPigment;
        	float closestAmb;
        	vec3 closestNormal;
        	float closestDiff;
        	float closestSpec;
        	float closestRoughness;
        	string type;
        	float t = -1; 
        	for(int i = 0; i < (signed) objects.size(); i++)
        	{
        		if(objects.at(i).type.compare("Sphere") == 0)
        		{
        			vec3 location = vec3(objects.at(i).location.x, objects.at(i).location.y, objects.at(i).location.z);
        			vec3 pigment = objects.at(i).pigment;
        			double radius = objects.at(i).length;
        			float ambient = objects.at(i).ambient;
        			float diffuse = objects.at(i).diffuse;
        			float specular = objects.at(i).specular;
        			float roughness = objects.at(i).roughness;
        			float a = dot(direction, direction);
        			float b = dot(2.f * direction , camera.location - location);
        			float c = dot(camera.location - location, camera.location - location) - (radius * radius);
        			float tSphere = quadraticSolver(a,b,c);
     
        			if(t == -1 && tSphere > 0)
        			{
        				t = tSphere;
        				closestPigment = pigment;
        				closestAmb = ambient;
        				closestDiff = diffuse;
        				vec3 point = camera.location + (tSphere * direction);
        				closestNormal = point - location;
        				closestNormal = normalize(closestNormal);
        				closestSpec = specular;
        				closestRoughness = roughness;
        				type = "Sphere";
        			}
        			else
        			{
        				if(tSphere < t && tSphere > 0 )
        				{

        					t = tSphere;
        					closestPigment = pigment;
        					type = "Sphere";
        					closestAmb = ambient;
        					closestDiff = diffuse;
        					vec3 point = camera.location + (tSphere * direction);
        					closestSpec = specular;
        					closestNormal = point - location;
        					closestNormal = normalize(closestNormal);
        					closestRoughness = roughness;
        				}
        			}
        		}
        		else if(objects.at(i).type.compare("Plane") == 0)
        		{
        			vec3 location = objects.at(i).location;
        			vec3 pigment = objects.at(i).pigment;
        			double length = objects.at(i).length;
        			float ambient = objects.at(i).ambient;
        			float diffuse = objects.at(i).diffuse;
        			float specular = objects.at(i).specular;
        			float roughness = objects.at(i).roughness;
        			float tPlane = (length - dot(camera.location, location))/dot(direction, location);
        			if(t == -1 && tPlane > 0)
        			{
        				t = tPlane;
        				closestPigment = pigment;
        				type = "Plane";
        				closestAmb = ambient;
        				closestDiff = diffuse;
        				closestNormal = location;
        				closestSpec = specular;
        				closestRoughness = roughness;
        			}
        			else
        			{
        					if(tPlane < t && tPlane > 0)
        					{
        						t = tPlane;
        						closestPigment = pigment;
        						type = "Plane";
        						closestAmb = ambient;
        						closestDiff = diffuse;	
        						closestNormal = location;
        						closestRoughness = roughness;
        						closestSpec = specular;
        					}
        			}
        		}
        	}
        cout << setprecision(4);
    	cout << "Pixel: [" << x << ", " << y << "] Ray: {" << camera.location.x << " " << camera.location.y << " "
    	<< camera.location.z << "} -> {" << direction.x << " " << direction.y << " " << direction.z << "}" << endl;
        	if(t == -1)
        	{
        		cout << "No Hit" << endl;
        	}
        	else
        	{

        	}

}*/

