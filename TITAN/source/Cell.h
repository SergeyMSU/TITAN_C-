#pragma once
#include"Header.h"


class Cell
{
public:
	vector<Yzel*> yzels;
	vector<Gran*> grans;
	int number = 0;                // номера начинаютс€ с единицы
	bool is_inner = false;         // явл€етс€ ли €чейка внутренней (которые считаютс€ отдельно)

	std::array < unordered_map<string, double>, 2> parameters;   // ѕараметры в €чейке (тут могут быть значени€
	// плазменных полей, значени€ дивергенций и т.д.
	// "rho", 'p', 'Vx', 'Vy', 'Vz', 'Bx', 'By', 'Bz'
	// 'Q' - ма€чок дл€ переноса и определени€ HP

	Cell* Get_Sosed(Gran* gr); // ѕолучить соседа через данную грань
	// возвращает nullptr, если соседа через данную грань нет, т.е. грань гранична€
	// ƒл€ работы функции номера €чеек должны быть актуальны

	friend Cell* Get_Sosed(Cell* C, Gran* gr);
	//  ак предыдуща€, но сравнение идЄт по указателю, не требует номеров €чеек

	double center[2][3];           // ÷ентр грани (также в предыдущий и следуюoий момент времени)
	double volume[2];

	void Culc_center(unsigned short int st_time);
	void Culc_volume(unsigned short int st_time, unsigned short int method = 1);

	double func_R(unsigned short int i_time); // –ассто€ние от центра до начала координат

};

