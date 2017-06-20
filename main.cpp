
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include "RayTrace.hpp"

using namespace std; 



int main(int argc, char* argv[])
{
	string command;
	string filename;
	
	shared_ptr<RayTrace> rTrace;
	rTrace = make_shared<RayTrace>();
	if(argc < 3)
	{
		cout << "Invalid number of arguments" << endl;
		return -1;
	}
	command = argv[1];
	filename = argv[2];
	if(command.compare("render") == 0)
	{
		int width = stoi(argv[3]);
		int height = stoi(argv[4]);
		rTrace->parseFile(filename);
		int ss = 1;
		int argument = 1;
		string type;
		if(argc == 6)
		{

			type = argv[5];
			
			if(type.substr(0, 4).compare("-ss=") == 0)
			{
				ss = stoi(type.substr(4, type.length()));
			}
			if(type.compare("-fresnel") == 0)
			{
				argument = 3;
			}
			if(type.compare("-sds") == 0)
			{
				argument = 0;
			}
			if(type.compare("-gi") == 0)
			{
				argument = 4;
			}
			if(type.compare("-caustic") == 0)
			{
				argument = 9;
			}
			if(type.compare("-sdscaustic") == 0)
			{
				argument = 10;
			}
			

	
			rTrace->drawScene(width, height, type, ss, argument);
		}
		else
		{

			rTrace->drawScene(width, height, "Blinn-Phong", ss, argument);
		}

	}
	else if(command.compare("pixeltrace") == 0)
	{
		rTrace->parseFile(filename);
		int width = stoi(argv[3]);
		int height = stoi(argv[4]);
		int x = stoi(argv[5]);
		int y = stoi(argv[6]);
		rTrace->pixelTrace(width, height, x, y);
	}
	else if(command.compare("sceneinfo") == 0)
	{
		rTrace->parseFile(filename);
		rTrace->printSceneInfo();
	}
	else if(command.compare("pixelray") == 0)
	{
		rTrace->parseFile(filename);
		int width = stoi(argv[3]);
		int height = stoi(argv[4]);
		int x = stoi(argv[5]);
		int y = stoi(argv[6]);
		rTrace->printRayInfo(width, height, x, y);
	}
	else if(command.compare("printrays") == 0)
	{
		rTrace->parseFile(filename);
		int width = stoi(argv[3]);
		int height = stoi(argv[4]);
		int x = stoi(argv[5]);
		int y = stoi(argv[6]);
		rTrace->printRays(width, height, x, y);
	}
	else if(command.compare("pixelcolor") == 0)
	{
		/*rTrace->parseFile(filename);
		int width = stoi(argv[3]);
		int height = stoi(argv[4]);
		int x = stoi(argv[5]);
		int y = stoi(argv[6]);
		string type;
		if(argc == 8)
		{
			type = argv[7];
			//rTrace->printpixColorInfo(width, height, x, y, type);
		}
		else
		{
			//rTrace->printpixColorInfo(width, height, x, y, "Blinn-Phong");
		}*/
		
	}
	else if(command.compare("firsthit") == 0)
	{
		rTrace->parseFile(filename);
		int width = stoi(argv[3]);
		int height = stoi(argv[4]);
		int x = stoi(argv[5]);
		int y = stoi(argv[6]);
		rTrace->printFirstHit(width, height, x, y); 
	}
	else
	{
		cout << "Invalid command." << endl << "Available commands include : render sceneinfo pixelray" << endl;
		return -1;
	}
	
	

	



	
}





