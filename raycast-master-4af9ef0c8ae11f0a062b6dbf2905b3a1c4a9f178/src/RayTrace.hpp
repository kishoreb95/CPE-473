#ifndef RAYTRACE_HPP
#define RAYTRACE_HPP

#include <iostream>
#include <string>
#include <stdio.h>
#include <fstream>
#include <glm/glm.hpp>
#include <vector>
#include <algorithm>
#include <glm/gtc/matrix_transform.hpp>

using namespace glm;
using namespace std;


class RayTrace
{
public:

	RayTrace(){};
	~RayTrace(){};
	int parseFile(string filename);
	glm::vec3 parseVector(string str, double *filter = 0);
	int parseObject(const string str, int index, string type);
	int parseCamera(const string str, int index);
	int parseLightSource(const string str, int index);
	void printRayInfo(int width, int height, int x, int y);
	void printSceneInfo();
	void printFirstHit(int width, int height, int x, int y);
	void pixelTrace(int width, int height, int x, int y);
	vec3 getImageCoordinate(int x, int y, int width, int height, float subpixelM, float subpixelN, int ss);
	float fresnel_reflectance(float ior, vec3 rayOrigin, vec3 ray);

	void printRays(int width, int height, int x, int y);

struct object {
	vec3 location;
	vec3 location1;
	vec3 location2;
	vec3 pigment;
	double length;
	double ambient;
	double diffuse;
	double specular;
	double roughness;
	double reflection = 0;
	float refraction = 0;
	double filter;
	vec3 center;
	float ior;
	vec3 scale; // 0
	vec3 rotate; // 1
	vec3 translate; // 2
	string order;
	mat4 ModelMatrix;
	mat4 invModelMatrix;
	string type;
	vec3 minimum;
	vec3 maximum;
	vec3 min;
	int objectIndex;
	vec3 max;

};
struct intersection
{
	string type;
	vec3 ray;
	vec3 intersectionPoint;
	int objectIndex;
	string objectType;
	float T = -1;
	vec3 rayOrigin;
	vec3 normal;
	vec3 transformedRay;
	vec3 reflectionRay;
	vec3 transmission;
	vec3 shadowRay;
	vec3 ambientVector;
	vec3 transformedOrigin;
	vec3 diffuseVector;
	vec3 specularVector;
	float ambient = 0;
	float diffuse = 0;
	float specular = 0;
	float reflection = 0;
	float refraction = 0;
	float roughness = 0;
	float filter = 0;
	float ior = 1;
	vec3 pigment;
	bool inshadow = false;
	vec3 color;
	struct object *obj;
} ;
struct bvh_node {
	struct bvh_node* left = NULL;
	struct bvh_node* right = NULL;
	bool isLeaf = true;
	vector<object> obj;
	vec3 minimum;
	vec3 maximum;
};
	struct bvh_node* root;
	void initializeBVH_Node(struct bvh_node* bvh);
	float getTBox(vec3 rayOrigin, vec3 ray, vec3 minimum, vec3 maximum);
	void bvh_traversal(vec3 rayOrigin, vec3 ray, struct bvh_node* bvh, struct intersection* inter);
	void makeBVHTree(struct bvh_node* bvh, vector<RayTrace::object> objs, int axis);
	void createBox(glm::vec3 location, glm::vec3 location1, glm::vec3 pigment, double ambient,double diffuse, double specular, double roughness, double reflection, double refraction, double ior, string type, double filter, vec3 sc, vec3 rot, vec3 trans, string order);
	void calculateBBox(struct bvh_node* bvh, vector<struct object> objs)
	{
		if(objs.size() == 0)
		{

		}
		
		else
		{
			
			std::vector<RayTrace::object> sortedByXAxis;
			std::vector<RayTrace::object> sortedByYAxis;
			std::vector<RayTrace::object> sortedByZAxis;

			sortedByXAxis = sort_objects_on_axis(objs, 0);
			sortedByYAxis = sort_objects_on_axis(objs, 1);
			sortedByZAxis = sort_objects_on_axis(objs, 2);
			
			bvh->minimum = vec3(sortedByXAxis.at(0).min.x, sortedByYAxis.at(0).min.y, sortedByZAxis.at(0).min.z);

			sortedByXAxis = sort_objects_on_axis(objs, 3);
			sortedByYAxis = sort_objects_on_axis(objs, 4);
			sortedByZAxis = sort_objects_on_axis(objs, 5);

			bvh->maximum = vec3(sortedByXAxis.at(0).max.x, sortedByYAxis.at(0).max.y, sortedByZAxis.at(0).max.z);


		
		}
		
	}
	vector<RayTrace::object> left_half_arr(std::vector<RayTrace::object> objs)
	{
		vector<RayTrace::object> leftHalf;
		for(int i = 0; i < (signed) objs.size() / 2; i++)
		{
			leftHalf.push_back(objs.at(i));
		}
		return leftHalf;
	}

	vector<RayTrace::object> right_half_arr(std::vector<RayTrace::object> objs)
	{
		vector<RayTrace::object> rightHalf;
		for(int i = (signed) objs.size() / 2; i < (signed) objs.size(); i++)
		{
			rightHalf.push_back(objs.at(i));
		}
		return rightHalf;
	}
	vector<RayTrace::object> sort_objects_on_axis(std::vector<RayTrace::object> objs, int axis)
	{
		if(axis == 0)
		{
			std::sort(objs.begin(), objs.end(), compareFunctionxAxis);
		}
		if(axis == 1)
		{
			std::sort(objs.begin(), objs.end(), compareFunctionyAxis);
		}
		if(axis == 2)
		{
			std::sort(objs.begin(), objs.end(), compareFunctionzAxis);
		}
		if(axis == 3)
		{
			std::sort(objs.begin(), objs.end(), maximumX);
		}
		if(axis == 4)
		{
			std::sort(objs.begin(), objs.end(), maximumY);
		}
		if(axis == 5)
		{
			std::sort(objs.begin(), objs.end(), maximumZ);
		}
		if(axis == 6)
		{
			std::sort(objs.begin(), objs.end(), centerX);
		}
		if(axis == 7)
		{
			std::sort(objs.begin(), objs.end(), centerY);
		}
		if(axis == 8)
		{
			std::sort(objs.begin(), objs.end(), centerZ);
		}
		return objs;
	}
	void calculateTransformedBoundingBox(struct object *sp);
	vector<vec3> sortVertices(vector<vec3> vert, int axis)
	{
    if(axis == 0)
    {
        std::sort(vert.begin(), vert.end(), sortX);
    }
    if(axis == 1)
    {
        std::sort(vert.begin(), vert.end(), sortY);
    }
    if(axis == 2)
    {
        std::sort(vert.begin(), vert.end(), sortZ);
    }
    return vert;
    
	}
	bool static maximumX(RayTrace::object obj1, RayTrace::object obj2)
	{
		return (obj1.max.x > obj2.max.x);
	}
	bool static maximumY(RayTrace::object obj1, RayTrace::object obj2)
	{
		return (obj1.max.y > obj2.max.y);
	}
	bool static maximumZ(RayTrace::object obj1, RayTrace::object obj2)
	{
		return (obj1.max.z <  obj2.max.z);
	}

	bool static compareFunctionxAxis(RayTrace::object obj1, RayTrace::object obj2)
	{
		return(obj1.min.x < obj2.min.x);
	}

	bool static compareFunctionyAxis(RayTrace::object obj1, RayTrace::object obj2)
	{
		return (obj1.min.y < obj2.min.y);
	}

	bool static compareFunctionzAxis(RayTrace::object obj1, RayTrace::object obj2)
	{
		return (obj1.min.z > obj2.min.z);
	}

	bool static centerX(RayTrace::object obj1, RayTrace::object obj2)
	{
		return (obj1.center.x < obj2.center.x);
	}

	bool static centerY(RayTrace::object obj1, RayTrace::object obj2)
	{
		return (obj1.center.y < obj2.center.y);
	}

	bool static centerZ(RayTrace::object obj1, RayTrace::object obj2)
	{
		return (obj1.center.z < obj2.center.z);
	}
	bool static sortX(vec3 obj1, vec3 obj2)
{
    return (obj1.x < obj2.x);
}
bool static sortY(vec3 obj1, vec3 obj2)
{
    return (obj1.y < obj2.y);
}
bool static sortZ(vec3 obj1, vec3 obj2)
{
    return (obj1.z < obj2.z);
}




struct camera {
	vec3 location;
	vec3 up;
	vec3 right;
	vec3 look_at;
};

struct light_source{
	vec3 location;
	vec3 pigment;
};


	float checkForIntersection(vec3 point, vec3 lvector);

	void createObject(glm::vec3 location, glm::vec3 pigment, float length, double ambient, double diffuse,double specular, double roughness, double reflection, double refraction, double ior, string type, double filter, vec3 scale, vec3 rotate, vec3 translate, string order);
	vec3 recTrace(vec3 r, vec3 reflectionVector, int depth);
	void createCamera(glm::vec3 location, glm::vec3 up, glm::vec3 right, glm::vec3 look_at);
	void createLightSource(glm::vec3 location, glm::vec3 pigment);
	void printpixColorInfo(int width,int height,int  x, int y, string typ);
	void createTriangle(glm::vec3 x, glm::vec3 y, glm::vec3 z, glm::vec3 pigment, float radius, double ambient, double diffuse, double specular, double roughness, string type, double filter);
	void printSpheres();
	vec3 traceAndPrint(vec3 rayOrigin, vec3 ray, int depth, vec3 final_color);
	void getIntersection(vec3 rayOrigin, vec3 direction, struct intersection *inter, std::vector<object> objs);
	void printCamera();
	void printLightSources();
	float tTriangle(vec3 a, vec3 b, vec3 c, vec3 direction, vec3 rayOrigin);
	float betaTriangle(vec3 a, vec3 b, vec3 c, vec3 direction, vec3 rayOrigin);
	void drawScene(int width, int height, string type, int ss, int arg);
	float gammaTriangle(vec3 a, vec3 b, vec3 c, vec3 direction, vec3 rayOrigin);
	float quadraticSolver(float a, float b, float c);
	vec3 blinnPhong(vec3 rayOrigin, vec3 direction, struct intersection *inter, int type );
	vec3 traceRay(vec3 rayOrigin, vec3 ray, int depth, int type);
	vector<vec3> generate_hemisphere_sample_points(int numPoints);
	vec3 alignSampleVector(vec3 sample, vec3 normal);
	

private:
	std::vector<object> objects;
	int numObjects = 0;
	struct camera camera;
	std::vector<light_source> light_source;
	std::vector<intersection> intersections;
	void addIntersection(struct intersection inter)
	{
		intersections.push_back(inter);
	}
	void addObject(struct object obj)
	{
		obj.objectIndex = numObjects;
		objects.push_back(obj);
		numObjects++;
	}
	
	void editCamera(struct camera cam)
	{
		camera = cam;
	}

	void addLightSource(struct light_source ls)
	{
		light_source.push_back(ls);
	}





	
};






#endif 
