#pragma once
#include "Header.h"
#include <vector>

class Gran;

enum class Type_yzel{
	TS,    // 0    узел на TS
	HP, // 1       узел на HP
	BS, // 2       узел на BS
	Us, // 3       обычный узел (не идентифицирован), 
	//такие есть на границе между зонами (если там нет поверхности разрыва)
	// точно есть в местах где можно продолжть HP и BS
	Zone_1,  // 4   узел до TS
	Zone_2,  // 5   узел между TS и HP
	Zone_3,  // 6   узел между HP и BS
	Zone_4  // 7   узел после BS
};

class Yzel
{
public:
	double coord[2][3];   // [Временной слой координат, три координаты]
	//Luch* luch;           // Луч на котором расположен узел (каждый узел может лежать только на одном луче)
	int number = 0;                // номера начинаются с единицы
	vector<Gran*> grans;
	Type_yzel type = Type_yzel::Us;  // по умолчанию создаём обычный узел
	bool is_inner = false;

	double velocity[3];
	int num_velocity = 0;

	unordered_map<string, double> parameters;   // Параметры в узле (для интерполяции надо)

	Yzel();
	Yzel(const double& a, const double& b, const double& c);

	double func_R(int i_time); // Расстояние до начала координат
	double func_Ryz(int i_time); // Расстояние до начала координат

	friend double Yzel_distance(Yzel* A, Yzel* B, int time); // расстояние между двумя узлами в пространстве
	friend double Yzel_distance_x(Yzel* A, Yzel* B, int time); // расстояние между двумя узлами в пространстве без учёта х координаты
};
