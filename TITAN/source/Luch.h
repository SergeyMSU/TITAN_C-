#pragma once

#include"Header.h"

using namespace std;


class Luch
{
public:
	string type;
	static Geo_param* geo;
	map<string, double> parameters;    // Параметры луча (могут быть свои для каждого типа)
	// "the" - угол в сферической ск от оси х - ось набегающего потока - зенитный угол
	// "phi" - угол в сферической ск - азимутальный угол
	//Yzel* start;            
	//Yzel* end;              
	vector<Yzel*> Yzels;
	vector<Yzel*> Yzels_opor;
	//vector<Luch*> Luch_soseds;    // Лучи-соседи (по ним легко искать соседние точки)

	void dvigenie(int i_time);
};
