#include "Setka.h"
#include "Header.h"
#include <iostream>
#include <algorithm>
#include <fstream>


using namespace std;

Setka::Setka()
{
    this->All_name_luch.push_back(&this->A_Luch);
    this->All_name_luch.push_back(&this->B_Luch);
    this->All_name_luch.push_back(&this->C_Luch);
    this->All_name_luch.push_back(&this->D_Luch);
    this->All_name_luch.push_back(&this->E_Luch);
    this->All_name_luch.push_back(&this->H_Luch);
    this->All_name_luch.push_back(&this->G_Luch);
}

void Setka::Set_geo()
{
    this->geo.phi = 0.0;
    this->geo.tetta0 = const_pi / 6;
    this->geo.tetta1 = const_pi - const_pi / 6;
    this->geo.tetta2 = const_pi * 120.0 / 180.0;
    this->geo.N1 = 20;
    this->geo.N2 = 10;
    this->geo.N3 = 15; //10
    this->geo.N4 = 10;
    this->geo.N5 = 8;

    this->geo.M0 = 40;//14;
    this->geo.M1 = 15;//15;
    this->geo.M11 = 3;
    this->geo.M2 = 25;//10;
    this->geo.M3 = 20;
    this->geo.M4 = 15;
    this->geo.MF = 5;//4;

    this->geo.R0 = 1.0;
    this->geo.R1 = 3.0;
    this->geo.R2 = 6.0;
    this->geo.R3 = 9.0;
    this->geo.R4 = 12.0;
    this->geo.R5 = 20.0;
    this->geo.L6 = -10.0;
    this->geo.L7 = -20.0;

    this->geo.dd1 = 4.0;    // У этих величин реализован автоматический подбор
    this->geo.dd2 = 2.0;    // У этих величин реализован автоматический подбор

    this->geo.dd3 = 2.2;
    this->geo.dd4 = 2.2;

    this->geo.dd5 = 2.3;
    this->geo.dd6 = 3.0;

    this->geo.dd7 = 3.8;
    this->geo.dd8 = 5.3;
}

void Setka::Construct_initial()
{
    double x, y, z, the, r;

    // Заполняем узлы на А-лучах
    for (int i = 0; i < this->geo.N1; i++)  
    {
        the = this->geo.tetta0 + i * (const_pi/2.0 - this->geo.tetta0)/(this->geo.N1 - 1);
        auto A = new Luch();
        A->phi = this->geo.phi;
        A->the = the;
        A->type = "A";
        this->All_Luch.push_back(A);
        this->A_Luch.push_back(A);

        for (int j = 0; j < this->geo.M0; j++)
        {
            r = this->geo.R0 + j * (this->geo.R1 - this->geo.R0) / this->geo.M0;
            x = r * cos(the);
            y = r * sin(the) * cos(A->phi);
            z = r * sin(the) * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // Опорная точка на радиусе R1 
        r = this->geo.R1;
        x = r * cos(the);
        y = r * sin(the) * cos(A->phi);
        z = r * sin(the) * sin(A->phi);
        auto Yz1 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz1);
        A->Yzels_opor.push_back(Yz1);
        this->All_Yzel.push_back(Yz1);
        Yz1->luch = A;
        
        // Точки от R1 до R2 (TS)
        for (int j = 0; j < this->geo.M1; j++)
        {
            r = this->geo.R1 + (j + 1) * (this->geo.R2 - this->geo.R1) / (this->geo.M1 + 1);
            x = r * cos(the);
            y = r * sin(the) * cos(A->phi);
            z = r * sin(the) * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // Опорная точка на радиусе R2 (TS) 
        r = this->geo.R2;
        x = r * cos(the);
        y = r * sin(the) * cos(A->phi);
        z = r * sin(the) * sin(A->phi);
        auto Yz2 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz2);
        A->Yzels_opor.push_back(Yz2);
        this->All_Yzel.push_back(Yz2);
        Yz2->luch = A;


        // Точки от R2 (TS) до R3 (HP)
        for (int j = 0; j < this->geo.M2; j++)
        {
            r = this->geo.R2 + (j + 1) * (this->geo.R3 - this->geo.R2) / (this->geo.M2 + 1);
            x = r * cos(the);
            y = r * sin(the) * cos(A->phi);
            z = r * sin(the) * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // Опорная точка на радиусе R3 (HP) 
        r = this->geo.R3;
        x = r * cos(the);
        y = r * sin(the) * cos(A->phi);
        z = r * sin(the) * sin(A->phi);
        auto Yz3 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz3);
        A->Yzels_opor.push_back(Yz3);
        this->All_Yzel.push_back(Yz3);
        Yz3->luch = A;

        // Точки от R3 (HP) до R4 (BS)
        for (int j = 0; j < this->geo.M3; j++)
        {
            r = this->geo.R3 + (j + 1) * (this->geo.R4 - this->geo.R3) / (this->geo.M3 + 1);
            x = r * cos(the);
            y = r * sin(the) * cos(A->phi);
            z = r * sin(the) * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // Опорная точка на радиусе R4  
        r = this->geo.R4;
        x = r * cos(the);
        y = r * sin(the) * cos(A->phi);
        z = r * sin(the) * sin(A->phi);
        auto Yz4 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz4);
        A->Yzels_opor.push_back(Yz4);
        this->All_Yzel.push_back(Yz4);
        Yz4->luch = A;


        // Точки от R4 (BS) до R5
        for (int j = 0; j < this->geo.M4; j++)
        {
            r = this->geo.R4 + (j + 1) * (this->geo.R5 - this->geo.R4) / (this->geo.M4 + 1);
            x = r * cos(the);
            y = r * sin(the) * cos(A->phi);
            z = r * sin(the) * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // Опорная точка на радиусе R5  
        r = this->geo.R5;
        x = r * cos(the);
        y = r * sin(the) * cos(A->phi);
        z = r * sin(the) * sin(A->phi);
        auto Yz5 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz5);
        A->Yzels_opor.push_back(Yz5);
        this->All_Yzel.push_back(Yz5);
        Yz5->luch = A;

    }

    // Заполняем узлы на B-лучах
    for (int i = 0; i < this->geo.N2; i++)
    {
        the = const_pi / 2.0 + (i + 1) * (this->geo.tetta2 - const_pi / 2.0) / (this->geo.N2);
        auto A = new Luch();
        A->phi = this->geo.phi;
        A->the = the;
        A->type = "B";
        this->All_Luch.push_back(A);
        this->B_Luch.push_back(A);

        for (int j = 0; j < this->geo.M0; j++)
        {
            r = this->geo.R0 + j * (this->geo.R1 - this->geo.R0) / this->geo.M0;
            x = r * cos(the);
            y = r * sin(the) * cos(A->phi);
            z = r * sin(the) * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // Опорная точка на радиусе R1 
        r = this->geo.R1;
        x = r * cos(the);
        y = r * sin(the) * cos(A->phi);
        z = r * sin(the) * sin(A->phi);
        auto Yz1 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz1);
        A->Yzels_opor.push_back(Yz1);
        this->All_Yzel.push_back(Yz1);
        Yz1->luch = A;

        // Точки от R1 до R2 (TS)
        for (int j = 0; j < this->geo.M1; j++)
        {
            r = this->geo.R1 + (j + 1) * (this->geo.R2 - this->geo.R1) / (this->geo.M1 + 1);
            x = r * cos(the);
            y = r * sin(the) * cos(A->phi);
            z = r * sin(the) * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // Опорная точка на радиусе R2 (TS) 
        r = this->geo.R2;
        x = r * cos(the);
        y = r * sin(the) * cos(A->phi);
        z = r * sin(the) * sin(A->phi);
        auto Yz2 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz2);
        A->Yzels_opor.push_back(Yz2);
        this->All_Yzel.push_back(Yz2);
        Yz2->luch = A;

        // Точки от R2 (TS) до R2 ++
        for (int j = 0; j < this->geo.M11; j++)
        {
            r = this->geo.R2 + (j + 1) * (this->geo.R3 - this->geo.R2) / (this->geo.M2 + 1);
            x = r * cos(the);
            y = r * sin(the) * cos(A->phi);
            z = r * sin(the) * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // Делаем непрямые лучи
        double x0, y0, x1, y1, t1, t2, tt1, tt2;

        t1 = cos(the) * this->geo.dd3;
        tt1 = sin(the) * this->geo.dd3;
        t2 = 0.0;
        tt2 = 1.0 * this->geo.dd4;
        x0 = x;        // Это координаты с последней итерации цикла
        y0 = y;
        x1 = x;
        y1 = this->geo.R3;

        double a1, b1, c1, d1, a2, b2, c2, d2;

        a1 = x0;
        b1 = t1;
        c1 = -2 * t1 - t2 - 3.0 * x0 + 3.0 * x1;
        d1 = t1 + t2 + 2.0 * x0 - 2.0 * x1;

        a2 = y0;
        b2 = tt1;
        c2 = -2 * tt1 - tt2 - 3.0 * y0 + 3.0 * y1;
        d2 = tt1 + tt2 + 2.0 * y0 - 2.0 * y1;

        // Точки от R2 (TS) до R3 (HP)
        for (int j = 0; j < this->geo.M2 - this->geo.M11; j++)
        {
            double s = 1.0 * (j + 1) / (this->geo.M2 - this->geo.M11 + 1);
            x = a1 + b1 * s + c1 * s * s + d1 * s * s * s;
            y = a2 + b2 * s + c2 * s * s + d2 * s * s * s;
            z = y * sin(A->phi);
            y = y * cos(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // Опорная точка на радиусе R3 (HP) 
        r = this->geo.R3;
        x = x0;
        y = r * cos(A->phi);
        z = r * sin(A->phi);
        auto Yz3 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz3);
        A->Yzels_opor.push_back(Yz3);
        this->All_Yzel.push_back(Yz3);
        Yz3->luch = A;

        // Точки от R3 (HP) до R4 (BS)
        for (int j = 0; j < this->geo.M3; j++)
        {
            r = this->geo.R3 + (j + 1) * (this->geo.R4 - this->geo.R3) / (this->geo.M3 + 1);
            x = x0;
            y = r * cos(A->phi);
            z = r * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // Опорная точка на радиусе R4  
        r = this->geo.R4;
        x = x0;
        y = r * cos(A->phi);
        z = r * sin(A->phi);
        auto Yz4 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz4);
        A->Yzels_opor.push_back(Yz4);
        this->All_Yzel.push_back(Yz4);
        Yz4->luch = A;


        // Точки от R4 (BS) до R5
        for (int j = 0; j < this->geo.M4; j++)
        {
            r = this->geo.R4 + (j + 1) * (this->geo.R5 - this->geo.R4) / (this->geo.M4 + 1);
            x = x0;
            y = r * cos(A->phi);
            z = r * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // Опорная точка на радиусе R4  
        r = this->geo.R5;
        x = x0;
        y = r * cos(A->phi);
        z = r * sin(A->phi);
        auto Yz5 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz5);
        A->Yzels_opor.push_back(Yz5);
        this->All_Yzel.push_back(Yz5);
        Yz5->luch = A;

    }

    // Заполняем узлы на C-лучах
    for (int i = 0; i < this->geo.N3; i++)
    {
        the = this->geo.tetta2 + (i + 1) * (this->geo.tetta1 - this->geo.tetta2) / (this->geo.N3);
        auto A = new Luch();
        A->phi = this->geo.phi;
        A->the = the;
        A->type = "C";
        this->All_Luch.push_back(A);
        this->C_Luch.push_back(A);

        for (int j = 0; j < this->geo.M0; j++)
        {
            r = this->geo.R0 + j * (this->geo.R1 - this->geo.R0) / this->geo.M0;
            x = r * cos(the);
            y = r * sin(the) * cos(A->phi);
            z = r * sin(the) * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // Опорная точка на радиусе R1 
        r = this->geo.R1;
        x = r * cos(the);
        y = r * sin(the) * cos(A->phi);
        z = r * sin(the) * sin(A->phi);
        auto Yz1 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz1);
        A->Yzels_opor.push_back(Yz1);
        this->All_Yzel.push_back(Yz1);
        Yz1->luch = A;

        // Точки от R1 до R2 (TS)
        for (int j = 0; j < this->geo.M1; j++)
        {
            r = this->geo.R1 + (j + 1) * (this->geo.R2 - this->geo.R1) / (this->geo.M1 + 1);
            x = r * cos(the);
            y = r * sin(the) * cos(A->phi);
            z = r * sin(the) * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // Опорная точка на радиусе R2 (TS) 
        r = this->geo.R2;
        x = r * cos(the);
        y = r * sin(the) * cos(A->phi);
        z = r * sin(the) * sin(A->phi);
        auto Yz2 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz2);
        A->Yzels_opor.push_back(Yz2);
        this->All_Yzel.push_back(Yz2);
        Yz2->luch = A;

        // Точки от R2 (TS) до R2 ++
        for (int j = 0; j < this->geo.M11; j++)
        {
            r = this->geo.R2 + (j + 1) * (this->geo.R3 - this->geo.R2) / (this->geo.M2 + 1); // Здесь не правильное 
            // расстояние взято  R3 тут не при чём. Но здесь млжно оставить, для начального построения пойдёт.
            x = r * cos(the);
            y = r * sin(the) * cos(A->phi);
            z = r * sin(the) * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // Делаем непрямые лучи
        double x0, y0, x1, y1, t1, t2, tt1, tt2;

        if (i == 0)
        {
            // Вот это не понятно откуда 
            t1 = (4 * cos(the) - 6 * 1.0) / 10.0 * this->geo.dd5; // Отклоняем вектор к хвосту немного сильнее
            tt1 = sin(the) * this->geo.dd5;
        }
        else
        {
            t1 = (cos(the) - 1.0) / 2.0 * this->geo.dd5;
            tt1 = sin(the) * this->geo.dd5;
        }
        t2 = -1.0 * this->geo.dd6;
        tt2 = 0.0;
        x0 = x;        // Это координаты с последней итерации цикла
        y0 = y;
        x1 = this->geo.L6;
        y1 = y;

        double a1, b1, c1, d1, a2, b2, c2, d2;

        a1 = x0;
        b1 = t1;
        c1 = -2 * t1 - t2 - 3.0 * x0 + 3.0 * x1;
        d1 = t1 + t2 + 2.0 * x0 - 2.0 * x1;

        a2 = y0;
        b2 = tt1;
        c2 = -2 * tt1 - tt2 - 3.0 * y0 + 3.0 * y1;
        d2 = tt1 + tt2 + 2.0 * y0 - 2.0 * y1;

        // Точки от R2 (TS) до L6
        for (int j = 0; j < this->geo.N4; j++)
        {
            double s = 1.0 * (j + 1) / (this->geo.N4 + 1);
            x = a1 + b1 * s + c1 * s * s + d1 * s * s * s;
            y = a2 + b2 * s + c2 * s * s + d2 * s * s * s;
            z = y * sin(A->phi);
            y = y * cos(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // Точки от L6 до L7
        for (int j = 0; j < this->geo.N5; j++)
        {
            double s = 1.0 * (j + 1) / (this->geo.N4 + 1);
            x = this->geo.L6 + (j) * (this->geo.L7 - this->geo.L6)/(this->geo.N5 - 1);
            y = y;
            z = z;
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;

            if(j == 0) A->Yzels_opor.push_back(Yz);
            if(j == this->geo.N5 - 1) A->Yzels_opor.push_back(Yz);
        }
    }

    // Заполняем узлы на D-лучах
    for (int i = 0; i < this->geo.N4 + this->geo.N5 - 3; i++)
    {
        the = 0.0;
        auto A = new Luch();
        A->phi = this->geo.phi;
        A->the = the;
        A->type = "D";
        this->All_Luch.push_back(A);
        this->D_Luch.push_back(A);

        double x, y, z;

        auto yz = this->C_Luch[0]->Yzels[this->geo.M0 + 1 + this->geo.M1 + 1 + this->geo.M11 + 3 + i];

        x = yz->coord[0][0];
        y = yz->coord[0][1];
        z = yz->coord[0][2];
        y = sqrt(y * y + z * z);
        A->Yzels.push_back(yz);
        A->Yzels_opor.push_back(yz);

        // Делаем непрямые лучи
        double x0, y0, x1, y1, t1, t2, tt1, tt2;

        t1 = 0.0;
        tt1 = 1.0 * 2;
        t2 = 0.0;
        tt2 = 1.0 * 2;
        x0 = x;        // Это координаты с последней итерации цикла
        y0 = y;
        x1 = x;
        y1 = this->geo.R3;

        double a1, b1, c1, d1, a2, b2, c2, d2;

        a1 = x0;
        b1 = t1;
        c1 = -2 * t1 - t2 - 3.0 * x0 + 3.0 * x1;
        d1 = t1 + t2 + 2.0 * x0 - 2.0 * x1;

        a2 = y0;
        b2 = tt1;
        c2 = -2 * tt1 - tt2 - 3.0 * y0 + 3.0 * y1;
        d2 = tt1 + tt2 + 2.0 * y0 - 2.0 * y1;

        // Точки от R2 (TS) до R3 (HP)
        for (int j = 0; j < this->geo.M2 - this->geo.M11; j++)
        {
            double s = 1.0 * (j + 1) / (this->geo.M2 - this->geo.M11 + 1);
            //x = a1 + b1 * s + c1 * s * s + d1 * s * s * s;
            //y = a2 + b2 * s + c2 * s * s + d2 * s * s * s;
            x = x0 + (x1 - x0) * (j + 1) / (this->geo.M2 - this->geo.M11 + 1);
            y = y0 + (y1 - y0) * (j + 1) / (this->geo.M2 - this->geo.M11 + 1);
            z = y * sin(A->phi);
            y = y * cos(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // Опорная точка на радиусе R3 (HP) 
        r = this->geo.R3;
        x = x0;
        y = r * cos(A->phi);
        z = r * sin(A->phi);
        auto Yz3 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz3);
        A->Yzels_opor.push_back(Yz3);
        this->All_Yzel.push_back(Yz3);
        Yz3->luch = A;

        // Точки от R3 (HP) до R4 (BS)
        for (int j = 0; j < this->geo.M3; j++)
        {
            r = this->geo.R3 + (j + 1) * (this->geo.R4 - this->geo.R3) / (this->geo.M3 + 1);
            x = x0;
            y = r * cos(A->phi);
            z = r * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // Опорная точка на радиусе R4  
        r = this->geo.R4;
        x = x0;
        y = r * cos(A->phi);
        z = r * sin(A->phi);
        auto Yz4 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz4);
        A->Yzels_opor.push_back(Yz4);
        this->All_Yzel.push_back(Yz4);
        Yz4->luch = A;


        // Точки от R4 (BS) до R5
        for (int j = 0; j < this->geo.M4; j++)
        {
            r = this->geo.R4 + (j + 1) * (this->geo.R5 - this->geo.R4) / (this->geo.M4 + 1);
            x = x0;
            y = r * cos(A->phi);
            z = r * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // Опорная точка на радиусе R4  
        r = this->geo.R5;
        x = x0;
        y = r * cos(A->phi);
        z = r * sin(A->phi);
        auto Yz5 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz5);
        A->Yzels_opor.push_back(Yz5);
        this->All_Yzel.push_back(Yz5);
        Yz5->luch = A;


    }

    // Заполняем узлы на E-лучах
    for (int i = 0; i < 4; i++)
    {
        the = 0.0;
        auto A = new Luch();
        A->phi = this->geo.phi;
        A->the = the;
        A->type = "E";
        this->All_Luch.push_back(A);
        this->E_Luch.push_back(A);

        auto yz1 = this->B_Luch[size(this->B_Luch) - 1]->Yzels_opor[2];
        auto yz2 = this->D_Luch[0]->Yzels_opor[1];

        double x, y, r;

        x = yz1->coord[0][0] + (i + 1) * (yz2->coord[0][0] - yz1->coord[0][0]) / 5;
        y = yz1->coord[0][1] + (i + 1) * (yz2->coord[0][1] - yz1->coord[0][1]) / 5;
        z = yz1->coord[0][2] + (i + 1) * (yz2->coord[0][2] - yz1->coord[0][2]) / 5;

        r = sqrt(y * y + z * z);

        // Опорная точка на радиусе R3 (HP) 
        r = this->geo.R3;
        y = r * cos(A->phi);
        z = r * sin(A->phi);
        auto Yz3 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz3);
        A->Yzels_opor.push_back(Yz3);
        this->All_Yzel.push_back(Yz3);
        Yz3->luch = A;

        // Точки от R3 (HP) до R4 (BS)
        for (int j = 0; j < this->geo.M3; j++)
        {
            r = this->geo.R3 + (j + 1) * (this->geo.R4 - this->geo.R3) / (this->geo.M3 + 1);
            y = r * cos(A->phi);
            z = r * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // Опорная точка на радиусе R4  
        r = this->geo.R4;
        y = r * cos(A->phi);
        z = r * sin(A->phi);
        auto Yz4 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz4);
        A->Yzels_opor.push_back(Yz4);
        this->All_Yzel.push_back(Yz4);
        Yz4->luch = A;


        // Точки от R4 (BS) до R5
        for (int j = 0; j < this->geo.M4; j++)
        {
            r = this->geo.R4 + (j + 1) * (this->geo.R5 - this->geo.R4) / (this->geo.M4 + 1);
            y = r * cos(A->phi);
            z = r * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // Опорная точка на радиусе R4  
        r = this->geo.R5;
        y = r * cos(A->phi);
        z = r * sin(A->phi);
        auto Yz5 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz5);
        A->Yzels_opor.push_back(Yz5);
        this->All_Yzel.push_back(Yz5);
        Yz5->luch = A;


    }

    // Добавляем H-лучи
    for (int i = 0; i < this->geo.MF - 1; i++)
    {
        the = 0.0;
        auto A = new Luch();
        A->phi = this->geo.phi;
        A->the = the;
        A->type = "H";
        this->All_Luch.push_back(A);
        this->H_Luch.push_back(A);

        int nn = this->B_Luch.size();
        auto yz1 = this->B_Luch[nn - 1]->Yzels[this->geo.M0 + this->geo.M1 + this->geo.M11 + 2 + i];
        auto yz3 = this->B_Luch[nn - 2]->Yzels[this->geo.M0 + this->geo.M1 + this->geo.M11 + 2 + i];
        auto yz2 = this->D_Luch[0]->Yzels[1 + i];
        auto yz4 = this->D_Luch[1]->Yzels[1 + i];

        A->Yzels.push_back(yz1);
        A->Yzels_opor.push_back(yz3); //
        A->Yzels_opor.push_back(yz1);

        // Делаем непрямые лучи
        double x0, y0, x1, y1, t1, t2, tt1, tt2;

        t1 = (yz1->coord[0][0] - yz3->coord[0][0]) * this->geo.dd7;
        tt1 = (yz1->coord[0][1] - yz3->coord[0][1]) * this->geo.dd7;
        t2 = (yz4->coord[0][0] - yz2->coord[0][0]) * this->geo.dd8;
        tt2 = (yz4->coord[0][1] - yz2->coord[0][1]) * this->geo.dd8;
        x0 = yz1->coord[0][0];        
        y0 = yz1->coord[0][1];
        x1 = yz2->coord[0][0];       
        y1 = yz2->coord[0][1];

        double a1, b1, c1, d1, a2, b2, c2, d2;

        a1 = x0;
        b1 = t1;
        c1 = -2 * t1 - t2 - 3.0 * x0 + 3.0 * x1;
        d1 = t1 + t2 + 2.0 * x0 - 2.0 * x1;

        a2 = y0;
        b2 = tt1;
        c2 = -2 * tt1 - tt2 - 3.0 * y0 + 3.0 * y1;
        d2 = tt1 + tt2 + 2.0 * y0 - 2.0 * y1;


        for (int j = 0; j < 4; j++)
        {
            double s = 1.0 * (j + 1) / (5);
            x = a1 + b1 * s + c1 * s * s + d1 * s * s * s;
            y = a2 + b2 * s + c2 * s * s + d2 * s * s * s;
            z = y * sin(A->phi);
            y = y * cos(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        A->Yzels.push_back(yz2);
        A->Yzels_opor.push_back(yz4);
        A->Yzels_opor.push_back(yz2);
    }

    // Добавляем G-лучи
    for (int i = 0; i < 4; i++)
    {
        the = 0.0;
        auto A = new Luch();
        A->phi = this->geo.phi;
        A->the = the;
        A->type = "G";
        this->All_Luch.push_back(A);
        this->G_Luch.push_back(A);

        int k = this->H_Luch.size();

        auto yz1 = this->H_Luch[k - 1]->Yzels[i + 1];
        auto yz2 = this->E_Luch[i]->Yzels_opor[0];
        auto yz3 = this->H_Luch[k - 2]->Yzels[i + 1];

        auto yz4 = this->B_Luch[this->B_Luch.size() - 1]->Yzels[this->geo.M0 + 1 + geo.M1 + 1 + geo.M2];
        auto yz5 = this->B_Luch[this->B_Luch.size() - 1]->Yzels[this->geo.M0 + 1 + geo.M1 + 1 + geo.M2 - 1];

        double x, y, r;

        A->Yzels.push_back(yz1);
        A->Yzels_opor.push_back(yz3);
        A->Yzels_opor.push_back(yz1);

        // Делаем непрямые лучи
        double x0, y0, x1, y1, t1, t2, tt1, tt2;

        // Для подбора коэффициентов
        double l1 = sqrt(kv((yz1->coord[0][0] - yz3->coord[0][0])) +
            kv((yz1->coord[0][1] - yz3->coord[0][1])));
        double l3 = sqrt(kv((yz4->coord[0][0] - yz5->coord[0][0])) +
            kv((yz4->coord[0][1] - yz5->coord[0][1])));

        double a1, b1, c1, d1, a2, b2, c2, d2;
        bool goto_aa1;

    aa1:
        goto_aa1 = false;
        t1 = (yz1->coord[0][0] - yz3->coord[0][0]) * this->geo.dd1;
        tt1 = (yz1->coord[0][1] - yz3->coord[0][1]) * this->geo.dd1;
        t2 = 0.0;
        tt2 = 1.0 * this->geo.dd2;
        x0 = yz1->coord[0][0];        
        y0 = yz1->coord[0][1];
        x1 = yz2->coord[0][0];
        y1 = yz2->coord[0][1];

        

        a1 = x0;
        b1 = t1;
        c1 = -2 * t1 - t2 - 3.0 * x0 + 3.0 * x1;
        d1 = t1 + t2 + 2.0 * x0 - 2.0 * x1;

        a2 = y0;
        b2 = tt1;
        c2 = -2 * tt1 - tt2 - 3.0 * y0 + 3.0 * y1;
        d2 = tt1 + tt2 + 2.0 * y0 - 2.0 * y1;

        // Точки от R2 (TS) до R3 (HP)
        for (int j = 0; j < this->geo.M2 - this->geo.M11 - this->geo.MF + 1; j++)
        {
            double s = 1.0 * (j + 1.0) / (this->geo.M2 - this->geo.M11 - this->geo.MF + 1 + 1.0);
            x = a1 + b1 * s + c1 * s * s + d1 * s * s * s;
            y = a2 + b2 * s + c2 * s * s + d2 * s * s * s;

            // Подбор геометрических коэффициентов
            if (i == 0 && j == 0) 
            {
                double l2 = sqrt(kv(x - x0) + kv(y - y0));
                if (fabs(l2 - l1) / l1 * 100.0 > 2.0)
                {
                    if (l2 < l1)
                    {
                        this->geo.dd1 *= 1.01;
                    }
                    else
                    {
                        this->geo.dd1 *= 0.99;
                    }
                    //cout << "1 " << l1 << " " << l2 << " " << this->geo.dd1 << endl;
                    //system("pause");
                    goto_aa1 = true;
                    break;
                }

                double xx, yy, ss;
                ss = 1.0 * (this->geo.M2 - this->geo.M11 - this->geo.MF + 1) / (this->geo.M2 - this->geo.M11 - this->geo.MF + 1 + 1.0);
                xx = a1 + b1 * ss + c1 * ss * ss + d1 * ss * ss * ss;
                yy = a2 + b2 * ss + c2 * ss * ss + d2 * ss * ss * ss;
                double l4 = sqrt(kv(xx - x1) + kv(yy - y1));
                if (fabs(l4 - l3) / l1 * 100.0 > 2.0)
                {
                    if (l4 < l3)
                    {
                        this->geo.dd2 *= 1.01;
                    }
                    else
                    {
                        this->geo.dd2 *= 0.99;
                    }
                    //cout << "2 " << l3 << " " << l4 << " " << this->geo.dd2 << endl;
                    //system("pause");
                    goto_aa1 = true;
                    break;
                }

            }


            z = y * sin(A->phi);
            y = y * cos(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        if (goto_aa1 == true)
        {
            goto aa1;
        }


        A->Yzels.push_back(yz2);
        A->Yzels_opor.push_back(yz2);
    }




    // Заполняем ячейки 
    // /////////////////////////////////////////////////////////////////////////////////////////////

    // Сразу будем сохранять ячейки в массивы ячеек (для дальнейшего быстрого поиска соседей)


    // Задаём массив А-ячеек (они так названы, потому что ограничены А-лучами)
    this->A_Cells.resize(this->geo.N1 - 1);
    for (unsigned i = 0; i < this->geo.N1 - 1; i++)
    {
        this->A_Cells[i].resize(this->A_Luch[0]->Yzels.size() - 1);
    }
    for (unsigned i = 0; i < this->A_Cells.size(); i++)
    {
        for (unsigned j = 0; j < this->A_Cells[0].size(); j++)
        {
            this->A_Cells[i][j] = nullptr;
        }
    }


    // А - лучи
    int n1 = this->A_Luch.size();
    if (n1 - 1 != this->geo.N1 - 1) cout << "Problem 1204846374563655463908" << endl;

    for (int i = 0; i < n1 - 1; i++)
    {
        int n2 = this->A_Luch[i]->Yzels.size();
        for (int j = 0; j < n2 - 1; j++)
        {
            auto A = new Cell();
            A->Yzels.push_back(this->A_Luch[i]->Yzels[j]);
            A->Yzels.push_back(this->A_Luch[i]->Yzels[j + 1]);
            A->Yzels.push_back(this->A_Luch[i + 1]->Yzels[j + 1]);
            A->Yzels.push_back(this->A_Luch[i + 1]->Yzels[j]);
            this->All_Cell.push_back(A);
            this->A_Cells[i][j] = A;
        }
    }


    // Задаём массив B-ячеек (они так названы, потому что ограничены B-лучами)
    this->B_Cells.resize(this->B_Luch.size() - 1);
    for (unsigned i = 0; i < this->B_Cells.size(); i++)
    {
        B_Cells[i].resize(this->B_Luch[0]->Yzels.size() - 1);
    }
    for (unsigned i = 0; i < B_Cells.size(); i++)
    {
        for (unsigned j = 0; j < B_Cells[0].size(); j++)
        {
            B_Cells[i][j] = nullptr;
        }
    }

    // B - лучи

    n1 = this->B_Luch.size();

    for (int i = 0; i < n1 - 1; i++)
    {
        int n2 = this->B_Luch[i]->Yzels.size();
        for (int j = 0; j < n2 - 1; j++)
        {
            auto A = new Cell();
            A->Yzels.push_back(this->B_Luch[i]->Yzels[j]);
            A->Yzels.push_back(this->B_Luch[i]->Yzels[j + 1]);
            A->Yzels.push_back(this->B_Luch[i + 1]->Yzels[j + 1]);
            A->Yzels.push_back(this->B_Luch[i + 1]->Yzels[j]);
            this->All_Cell.push_back(A);
            this->B_Cells[i][j] = A;
        }
    }


    // Задаём массив C-ячеек (они так названы, потому что ограничены C-лучами)
    this->C_Cells.resize(this->C_Luch.size() - 1);
    for (unsigned i = 0; i < this->C_Cells.size(); i++)
    {
        C_Cells[i].resize(this->C_Luch[0]->Yzels.size() - 1);
    }
    for (unsigned i = 0; i < C_Cells.size(); i++)
    {
        for (unsigned j = 0; j < C_Cells[0].size(); j++)
        {
            C_Cells[i][j] = nullptr;
        }
    }


    // C - лучи

    n1 = this->C_Luch.size();

    for (int i = 0; i < n1 - 1; i++)
    {
        int n2 = this->C_Luch[i]->Yzels.size();
        for (int j = 0; j < n2 - 1; j++)
        {
            auto A = new Cell();
            A->Yzels.push_back(this->C_Luch[i]->Yzels[j]);
            A->Yzels.push_back(this->C_Luch[i]->Yzels[j + 1]);
            A->Yzels.push_back(this->C_Luch[i + 1]->Yzels[j + 1]);
            A->Yzels.push_back(this->C_Luch[i + 1]->Yzels[j]);
            this->All_Cell.push_back(A);
            this->C_Cells[i][j] = A;
        }
    }

    // Задаём массив D-ячеек (они так названы, потому что ограничены D-лучами)
    this->D_Cells.resize(this->D_Luch.size() - 1);
    for (unsigned i = 0; i < this->D_Cells.size(); i++)
    {
        D_Cells[i].resize(this->D_Luch[0]->Yzels.size() - 1);
    }
    for (unsigned i = 0; i < D_Cells.size(); i++)
    {
        for (unsigned j = 0; j < D_Cells[0].size(); j++)
        {
            D_Cells[i][j] = nullptr;
        }
    }


    // D - лучи

    n1 = this->D_Luch.size();

    for (int i = 0; i < n1 - 1; i++)
    {
        int n2 = this->D_Luch[i]->Yzels.size();
        for (int j = 0; j < n2 - 1; j++)
        {
            auto A = new Cell();
            A->Yzels.push_back(this->D_Luch[i]->Yzels[j]);
            A->Yzels.push_back(this->D_Luch[i]->Yzels[j + 1]);
            A->Yzels.push_back(this->D_Luch[i + 1]->Yzels[j + 1]);
            A->Yzels.push_back(this->D_Luch[i + 1]->Yzels[j]);
            this->All_Cell.push_back(A);
            this->D_Cells[i][j] = A;
        }
    }


    // Задаём массив E-ячеек (они так названы, потому что ограничены E-лучами)
    this->E_Cells.resize(this->E_Luch.size() - 1);
    for (unsigned i = 0; i < this->E_Cells.size(); i++)
    {
        E_Cells[i].resize(this->E_Luch[0]->Yzels.size() - 1);
    }
    for (unsigned i = 0; i < E_Cells.size(); i++)
    {
        for (unsigned j = 0; j < E_Cells[0].size(); j++)
        {
            E_Cells[i][j] = nullptr;
        }
    }

    // E - лучи

    n1 = this->E_Luch.size();

    for (int i = 0; i < n1 - 1; i++)
    {
        int n2 = this->E_Luch[i]->Yzels.size();
        for (int j = 0; j < n2 - 1; j++)
        {
            auto A = new Cell();
            A->Yzels.push_back(this->E_Luch[i]->Yzels[j]);
            A->Yzels.push_back(this->E_Luch[i]->Yzels[j + 1]);
            A->Yzels.push_back(this->E_Luch[i + 1]->Yzels[j + 1]);
            A->Yzels.push_back(this->E_Luch[i + 1]->Yzels[j]);
            this->All_Cell.push_back(A);
            this->E_Cells[i][j] = A;
        }
    }

    // Задаём массив G-ячеек (они так названы, потому что ограничены G-лучами)
    this->G_Cells.resize(this->G_Luch.size() - 1);
    for (unsigned i = 0; i < this->G_Cells.size(); i++)
    {
        G_Cells[i].resize(this->G_Luch[0]->Yzels.size() - 1);
    }
    for (unsigned i = 0; i < G_Cells.size(); i++)
    {
        for (unsigned j = 0; j < G_Cells[0].size(); j++)
        {
            G_Cells[i][j] = nullptr;
        }
    }


    // G - лучи

    n1 = this->G_Luch.size();

    for (int i = 0; i < n1 - 1; i++)
    {
        int n2 = this->G_Luch[i]->Yzels.size();
        for (int j = 0; j < n2 - 1; j++)
        {
            auto A = new Cell();
            A->Yzels.push_back(this->G_Luch[i]->Yzels[j]);
            A->Yzels.push_back(this->G_Luch[i]->Yzels[j + 1]);
            A->Yzels.push_back(this->G_Luch[i + 1]->Yzels[j + 1]);
            A->Yzels.push_back(this->G_Luch[i + 1]->Yzels[j]);
            this->All_Cell.push_back(A);
            this->G_Cells[i][j] = A;
        }
    }


    // Задаём массив H-ячеек (они так названы, потому что ограничены H-лучами)
    this->H_Cells.resize(this->H_Luch.size() - 1);
    for (unsigned i = 0; i < this->H_Cells.size(); i++)
    {
        H_Cells[i].resize(this->H_Luch[0]->Yzels.size() - 1);
    }
    for (unsigned i = 0; i < H_Cells.size(); i++)
    {
        for (unsigned j = 0; j < H_Cells[0].size(); j++)
        {
            H_Cells[i][j] = nullptr;
        }
    }


    // H - лучи

    n1 = this->H_Luch.size();

    for (int i = 0; i < n1 - 1; i++)
    {
        int n2 = this->H_Luch[i]->Yzels.size();
        for (int j = 0; j < n2 - 1; j++)
        {
            auto A = new Cell();
            A->Yzels.push_back(this->H_Luch[i]->Yzels[j]);
            A->Yzels.push_back(this->H_Luch[i + 1]->Yzels[j]);
            A->Yzels.push_back(this->H_Luch[i + 1]->Yzels[j + 1]);
            A->Yzels.push_back(this->H_Luch[i]->Yzels[j + 1]);
            this->All_Cell.push_back(A);
            this->H_Cells[i][j] = A;
        }
    }


    // BC - склейка
    for (int j = 0; j < this->geo.M0 + this->geo.M1 + this->geo.M11 + 1; j++)
    {
        int k = this->B_Luch.size();
        auto A = new Cell();
        A->Yzels.push_back(this->B_Luch[k - 1]->Yzels[j]);
        A->Yzels.push_back(this->B_Luch[k - 1]->Yzels[j + 1]);
        A->Yzels.push_back(this->C_Luch[0]->Yzels[j + 1]);
        A->Yzels.push_back(this->C_Luch[0]->Yzels[j]);
        this->All_Cell.push_back(A);
    }

    // AB - склейка
    for (int j = 0; j < this->B_Luch[0]->Yzels.size() - 1; j++)
    {
        int k = this->A_Luch.size();
        auto A = new Cell();
        A->Yzels.push_back(this->A_Luch[k - 1]->Yzels[j]);
        A->Yzels.push_back(this->A_Luch[k - 1]->Yzels[j + 1]);
        A->Yzels.push_back(this->B_Luch[0]->Yzels[j + 1]);
        A->Yzels.push_back(this->B_Luch[0]->Yzels[j]);
        this->All_Cell.push_back(A);
    }

    // BE - склейка
    for (int j = 0; j < this->E_Luch[0]->Yzels.size() - 1; j++)
    {
        int k = this->B_Luch.size();
        int nn = this->geo.M0 + this->geo.M1 + this->geo.M2 + 2;
        auto A = new Cell();
        A->Yzels.push_back(this->B_Luch[k - 1]->Yzels[nn + j]);
        A->Yzels.push_back(this->B_Luch[k - 1]->Yzels[nn + j + 1]);
        A->Yzels.push_back(this->E_Luch[0]->Yzels[j + 1]);
        A->Yzels.push_back(this->E_Luch[0]->Yzels[j]);
        this->All_Cell.push_back(A);
    }

    // ED - склейка
    for (int j = 0; j < this->E_Luch[0]->Yzels.size() - 1; j++)
    {
        int k = this->E_Luch.size();
        int nn = this->geo.M2 - this->geo.M11 + 1;
        auto A = new Cell();
        A->Yzels.push_back(this->E_Luch[k - 1]->Yzels[j]);
        A->Yzels.push_back(this->E_Luch[k - 1]->Yzels[j + 1]);
        A->Yzels.push_back(this->D_Luch[0]->Yzels[nn + j + 1]);
        A->Yzels.push_back(this->D_Luch[0]->Yzels[nn + j]);
        this->All_Cell.push_back(A);
    }


    // BG - склейка
    for (int j = 0; j < this->G_Luch[0]->Yzels.size() - 1; j++)
    {
        int k = this->B_Luch.size();
        int nn = this->geo.M0 + this->geo.M1 + this->geo.M11 + this->geo.MF;
        auto A = new Cell();
        A->Yzels.push_back(this->B_Luch[k - 1]->Yzels[nn + j]);
        A->Yzels.push_back(this->B_Luch[k - 1]->Yzels[nn + j + 1]);
        A->Yzels.push_back(this->G_Luch[0]->Yzels[j + 1]);
        A->Yzels.push_back(this->G_Luch[0]->Yzels[j]);
        this->All_Cell.push_back(A);
    }

    // GD - склейка
    for (int j = 0; j < this->G_Luch[0]->Yzels.size() - 1; j++)
    {
        int k = this->G_Luch.size();
        int nn = this->geo.MF - 1;
        auto A = new Cell();
        A->Yzels.push_back(this->G_Luch[k - 1]->Yzels[j]);
        A->Yzels.push_back(this->G_Luch[k - 1]->Yzels[j + 1]);
        A->Yzels.push_back(this->D_Luch[0]->Yzels[nn + j + 1]);
        A->Yzels.push_back(this->D_Luch[0]->Yzels[nn + j]);
        this->All_Cell.push_back(A);
    }

    // Ручная склейка особых точек

    if (true)
    {
        //1
        int k = this->B_Luch.size();
        auto yz1 = this->B_Luch[k - 1]->Yzels[this->geo.M0 + this->geo.M1 + this->geo.M11 + 1];
        auto yz2 = this->B_Luch[k - 1]->Yzels[this->geo.M0 + this->geo.M1 + this->geo.M11 + 2];
        auto yz3 = this->H_Luch[0]->Yzels[1];
        auto yz4 = this->C_Luch[0]->Yzels[this->geo.M0 + this->geo.M1 + this->geo.M11 + 1];

        auto A = new Cell();
        A->Yzels.push_back(yz1);
        A->Yzels.push_back(yz2);
        A->Yzels.push_back(yz3);
        A->Yzels.push_back(yz4);
        this->All_Cell.push_back(A);



        //3
        for (int i = 0; i < 3; i++)
        {
            yz1 = this->C_Luch[0]->Yzels[this->geo.M0 + this->geo.M1 + this->geo.M11 + 1 + i];
            yz2 = this->H_Luch[0]->Yzels[1 + i];
            yz3 = this->H_Luch[0]->Yzels[2 + i];
            yz4 = this->C_Luch[0]->Yzels[this->geo.M0 + this->geo.M1 + this->geo.M11 + 2 + i];

            A = new Cell();
            A->Yzels.push_back(yz1);
            A->Yzels.push_back(yz2);
            A->Yzels.push_back(yz3);
            A->Yzels.push_back(yz4);
            this->All_Cell.push_back(A);
        }

        //4 
        yz1 = this->C_Luch[0]->Yzels[this->geo.M0 + this->geo.M1 + this->geo.M11 + 1 + 3];
        yz2 = this->H_Luch[0]->Yzels[4];
        yz3 = this->D_Luch[0]->Yzels[1];
        yz4 = this->C_Luch[0]->Yzels[this->geo.M0 + this->geo.M1 + this->geo.M11 + 2 + 3];

        A = new Cell();
        A->Yzels.push_back(yz1);
        A->Yzels.push_back(yz2);
        A->Yzels.push_back(yz3);
        A->Yzels.push_back(yz4);
        this->All_Cell.push_back(A);
    }

    // Заполняем соседей ячеек
    // Сначала пробегаемся внутри массивов ячеек, потом будем между массивами склеивать алгоритмом поиска

    vector<vector<Cell*>>& Cells_ = this->A_Cells;
    int nn1, mm1;

    vector<vector<vector<Cell*>>*> all_Cells_;
    all_Cells_.push_back(&this->A_Cells);
    all_Cells_.push_back(&this->B_Cells);
    all_Cells_.push_back(&this->C_Cells);
    all_Cells_.push_back(&this->D_Cells);
    all_Cells_.push_back(&this->E_Cells);
    all_Cells_.push_back(&this->H_Cells);
    all_Cells_.push_back(&this->G_Cells);

    for (auto& ssCells_ : all_Cells_)
    {
        Cells_ = *ssCells_;
        nn1 = Cells_.size();
        mm1 = Cells_[0].size();
        for (int i = 0; i < nn1; i++)
        {
            for (int j = 0; j < mm1; j++)
            {
                int ip, im, jp, jm;
                ip = i + 1;
                im = i - 1;
                jp = j + 1;
                jm = j - 1;
                if (ip < nn1) Cells_[i][j]->Soseds.push_back(Cells_[ip][j]);
                if (im >= 0) Cells_[i][j]->Soseds.push_back(Cells_[im][j]);
                if (jp < mm1) Cells_[i][j]->Soseds.push_back(Cells_[i][jp]);
                if (jm >= 0) Cells_[i][j]->Soseds.push_back(Cells_[i][jm]);
            }
        }
    }

    // Перед тем, как склеивать, пробежимся по ячейкам и добавим их в узел
    for (auto& i : this->All_Cell)
    {
        for (auto& j : i->Yzels)
        {
            j->Cells.push_back(i);
        }
    }

    // Теперь будем склеивать ячейки через узлы
    for (auto& i : this->All_Yzel)
    {
        for (auto& j : i->Cells)
        {
            if (j->Soseds.size() == 4) continue;
            for (auto& jj : i->Cells)
            {
                if (jj == j) continue;

                // Проверям, может быть ячейки уже являются найденными соседями
                if (find(j->Soseds.begin(), j->Soseds.end(), jj) != j->Soseds.end()) continue;

                if (Cell_is_soseds(j, jj))
                {
                    j->Soseds.push_back(jj);
                    jj->Soseds.push_back(j);
                }
            }
        }
    }

}

void Setka::All_numerate()
{
    int k = 1;
    for (auto& i : this->All_Cell)
    {
        i->number = k;
        k++;
    }

    k = 1;
    for (auto& i : this->All_Yzel)
    {
        i->number = k;
        k++;
    }
}

bool Setka::Cell_is_soseds(Cell* A, Cell* B)
{
    int k = 0;
    for (auto& jj : B->Yzels)
    {
        if (find(A->Yzels.begin(), A->Yzels.end(), jj) != A->Yzels.end()) k++;
    }

    if (k == 2) return true;
    return false;
}

void Setka::Print_yzel()
{
    ofstream fout;
    fout.open("All_yzel.txt");

    for (auto& i : this->All_Yzel)
    {
        fout << i->coord[0][0] << " " << i->coord[0][1] << " " << i->coord[0][2] << endl;
    }
    fout.close();
}

void Setka::Print_yzels_opor(string name)
{
    ofstream fout;
    fout.open(name + "_All_opor_yzel.txt");

    vector<Luch*>* Luch_all = nullptr;

    if (name == "A") Luch_all = &this->A_Luch;
    if (name == "B") Luch_all = &this->B_Luch;
    if (name == "C") Luch_all = &this->C_Luch;
    if (name == "D") Luch_all = &this->D_Luch;
    if (name == "E") Luch_all = &this->E_Luch;
    if (name == "H") Luch_all = &this->H_Luch;
    if (name == "G") Luch_all = &this->G_Luch;

    for (auto& i : *Luch_all)
    {
        for (auto& j : i->Yzels_opor)
        {
            fout << j->coord[0][0] << " " << j->coord[0][1] << endl;
        }
    }

    fout.close();
}

void Setka::Print_cell()
{
    int ll = this->All_Cell.size();
    ofstream fout;
    fout.open("Setka_all_cell.txt");
    fout << "TITLE = \"HP\" ";
    fout << " VARIABLES = \"X\", \"Y\"  ZONE T= \"HP\", N=  " << 4 * ll;
    fout << " , E= " << 4 * ll;
    fout << " , F=FEPOINT, ET=LINESEG  " << endl;
    for (auto& i : this->All_Cell)
    {
        for (auto& j : i->Yzels)
        {
            fout << j->coord[0][0] << " " << j->coord[0][1] << endl;
        }
    }
    for (int i = 0; i < ll; i++)
    {
        fout << 4 * i + 1 << " " << 4 * i + 2 << endl;
        fout << 4 * i + 2 << " " << 4 * i + 3 << endl;
        fout << 4 * i + 3 << " " << 4 * i + 4 << endl;
        fout << 4 * i + 4 << " " << 4 * i + 1 << endl;
    }

    fout.close();
}

void Setka::Print_luch()
{
    int ll = 0;
    int nn = 0;
    for (auto& i : this->A_Luch)
    {
        ll += (i->Yzels.size() - 1);
        nn += i->Yzels.size();
    }
    for (auto& i : this->B_Luch)
    {
        ll += (i->Yzels.size() - 1);
        nn += i->Yzels.size();
    }
    for (auto& i : this->C_Luch)
    {
        ll += (i->Yzels.size() - 1);
        nn += i->Yzels.size();
    }
    for (auto& i : this->D_Luch)
    {
        ll += (i->Yzels.size() - 1);
        nn += i->Yzels.size();
    }
    for (auto& i : this->E_Luch)
    {
        ll += (i->Yzels.size() - 1);
        nn += i->Yzels.size();
    }
    for (auto& i : this->H_Luch)
    {
        ll += (i->Yzels.size() - 1);
        nn += i->Yzels.size();
    }
    for (auto& i : this->G_Luch)
    {
        ll += (i->Yzels.size() - 1);
        nn += i->Yzels.size();
    }


    ofstream fout;
    fout.open("Setka_all_luch.txt");
    fout << "TITLE = \"HP\" ";
    fout << " VARIABLES = \"X\", \"Y\", \"Num\"  ZONE T= \"HP\", N=  " << nn;
    fout << " , E= " << ll;
    fout << " , F=FEPOINT, ET=LINESEG  " << endl;


    for (auto& i : this->All_Luch)
    {
        int tt = 0;
        if (i->type == "A") tt = 1;
        if (i->type == "B") tt = 2;
        if (i->type == "C") tt = 3;
        if (i->type == "D") tt = 4;
        if (i->type == "E") tt = 5;
        if (i->type == "H") tt = 6;
        if (i->type == "G") tt = 7;
        for (auto& j : i->Yzels)
        {
            fout << j->coord[0][0] << " " << j->coord[0][1] << " " << tt << endl;
        }
    }

    int kk = 0;

    for (auto& i : this->All_Luch)
    {
        for (int j = 0; j < i->Yzels.size() - 1; j++)
        {
            fout << kk + j + 1 << " " << kk + j + 2 << endl;
        }
        kk += i->Yzels.size();
    }


    fout.close();
}

void Setka::Print_cell_soseds()
{
    int ll = 0;
    int nn = 0;
    for (auto& i : this->All_Cell)
    {
        ll += i->Soseds.size();
        i->Set_center();
    }
    nn = 2 * ll;

    ofstream fout;
    fout.open("Setka_all_soseds.txt");
    fout << "TITLE = \"HP\" ";
    fout << " VARIABLES = \"X\", \"Y\"  ZONE T= \"HP\", N=  " << nn;
    fout << " , E= " << ll;
    fout << " , F=FEPOINT, ET=LINESEG  " << endl;

    double x0, y0, x1, y1;

    for (auto& i : this->All_Cell)
    {
        x0 = i->xc;
        y0 = i->yc;

        for (auto& j : i->Soseds)
        {
            x1 = j->xc;
            y1 = j->yc;
            fout << x0 << " " << y0 << endl;
            fout << (x0 + x1)/2.0 << " " << (y0 + y1) / 2.0 << endl;
        }
    }

    for (int i = 0; i < ll; i++)
    {
        fout << 2 * i + 1 << " " << 2 * i + 2 << endl;
    }


    fout.close();
}

void Setka::Print_cell_center()
{
    for (auto& i : this->All_Cell)
    {
        i->Set_center();
    }

    ofstream fout;
    fout.open("All_Cell_center.txt");

   for (auto& i : this->All_Cell)
    {
        fout << i->xc << " " << i->yc << endl;
    }
    fout.close();
}

void Setka::Print_cell_center_test()
{
    for (auto& i : this->All_Cell)
    {
        i->Set_center();
    }

    ofstream fout;
    fout.open("All_Cell_center_test.txt");

    for (auto& i : this->All_Cell)
    {
        if (i->Soseds.size() == 4) continue;
        if (i->Soseds.size() <= 1) cout << "TEST ERROR  09328764t97364269 = " << i->Soseds.size() <<  endl;
        if (i->Soseds.size() > 4) cout << "TEST ERROR  07896856765874645676755 = " << i->Soseds.size() <<  endl;
        
        //if (i->Soseds.size() != 3) continue;
        //fout << i->xc << " " << i->yc << endl;
    }
    fout.close();
}

void Setka::Save_for_3D(string name)
{
    this->All_numerate();

    // Пишем информационный файл для сетки (его не надо будет считывать программой, он для пользователя)
    ofstream fout;
    fout.open(name + "_2D_Setka_INFO.txt");
    fout << "Vsego yacheek v setke = " << this->All_Cell.size() << endl;
    fout << "Vsego uzlov v setke = " << this->All_Yzel.size() << endl;
    fout << "Vsego luchej v setke = " << this->All_Luch.size() << endl;
    fout << "Luchej tipa A v setke = " << this->A_Luch.size() << endl;
    fout << "Luchej tipa B v setke = " << this->B_Luch.size() << endl;
    fout << "Luchej tipa C v setke = " << this->C_Luch.size() << endl;
    fout << "Luchej tipa D v setke = " << this->D_Luch.size() << endl;
    fout << "Luchej tipa E v setke = " << this->E_Luch.size() << endl;
    fout << "Luchej tipa H v setke = " << this->H_Luch.size() << endl;
    fout << "Luchej tipa G v setke = " << this->G_Luch.size() << endl;
    fout << "Geometriya ----------------------------------------------------" << endl;
    fout << "tetta0 = " << this->geo.tetta0 << endl;
    fout << "tetta1 = " << this->geo.tetta1 << endl;
    fout << "tetta2 = " << this->geo.tetta2 << endl;
    fout << "N1 (CHislo A - napravlyayushchih) = " << this->geo.N1 << endl;
    fout << "N2 (CHislo B - napravlyayushchih) = " << this->geo.N2 << endl;
    fout << "N3 (CHislo C - napravlyayushchih) = " << this->geo.N3 << endl;
    fout << "N4 (CHislo D - napravlyayushchih v podvizhnoj zone + 3 (sdvig iz-za drugih napravlyayushchih)) = " << this->geo.N4 << endl;
    fout << "N5 CHislo D - napravlyayushchih (v nepodvizhnoj zone, levee L6) = " << this->geo.N5 << endl;
    fout << "M0 = " << this->geo.M0 << endl;
    fout << "M1 = " << this->geo.M1 << endl;
    fout << "M11 = " << this->geo.M11 << endl;
    fout << "M2 = " << this->geo.M2 << endl;
    fout << "M3 = " << this->geo.M3 << endl;
    fout << "M4 = " << this->geo.M4 << endl;
    fout << "MF = " << this->geo.MF << endl;
    fout << "R0 = " << this->geo.R0 << endl;
    fout << "R1 = " << this->geo.R1 << endl;
    fout << "R2 = " << this->geo.R2 << endl;
    fout << "R3 = " << this->geo.R3 << endl;
    fout << "R4 = " << this->geo.R4 << endl;
    fout << "R5 = " << this->geo.R5 << endl;
    fout << "L6 = " << this->geo.L6 << endl;
    fout << "L7 = " << this->geo.L7 << endl;
    fout << "dd1 = " << this->geo.dd1 << endl;
    fout << "dd2 = " << this->geo.dd2 << endl;
    fout << "dd3 = " << this->geo.dd3 << endl;
    fout << "dd4 = " << this->geo.dd4 << endl;
    fout << "dd5 = " << this->geo.dd5 << endl;
    fout << "dd6 = " << this->geo.dd6 << endl;
    fout << "dd7 = " << this->geo.dd7 << endl;
    fout << "dd8 = " << this->geo.dd8 << endl;
    fout.close();

    ofstream out(name + "_2D_Setka.bin", ios::binary | ios::out);
    int n;
    double x, y;

    // Записываем номер версии файла вывода (для того, чтобы можно было отслеживать версии)
    n = 2;
    out.write((char*)&n, sizeof n);

    // Записываем параметры сетки   this->geo
    if (true)
    {
        x = this->geo.tetta0;
        out.write((char*)&x, sizeof x);
        x = this->geo.tetta1;
        out.write((char*)&x, sizeof x);
        x = this->geo.tetta2;
        out.write((char*)&x, sizeof x);

        n = this->geo.N1;
        out.write((char*)&n, sizeof n);
        n = this->geo.N2;
        out.write((char*)&n, sizeof n);
        n = this->geo.N3;
        out.write((char*)&n, sizeof n);
        n = this->geo.N4;
        out.write((char*)&n, sizeof n);
        n = this->geo.N5;
        out.write((char*)&n, sizeof n);

        n = this->geo.M0;
        out.write((char*)&n, sizeof n);
        n = this->geo.M1;
        out.write((char*)&n, sizeof n);
        n = this->geo.M11;
        out.write((char*)&n, sizeof n);
        n = this->geo.M2;
        out.write((char*)&n, sizeof n);
        n = this->geo.M3;
        out.write((char*)&n, sizeof n);
        n = this->geo.M4;
        out.write((char*)&n, sizeof n);
        n = this->geo.MF;
        out.write((char*)&n, sizeof n);

        x = this->geo.R0;
        out.write((char*)&x, sizeof x);
        x = this->geo.R1;
        out.write((char*)&x, sizeof x);
        x = this->geo.R2;
        out.write((char*)&x, sizeof x);
        x = this->geo.R3;
        out.write((char*)&x, sizeof x);
        x = this->geo.R4;
        out.write((char*)&x, sizeof x);
        x = this->geo.R5;
        out.write((char*)&x, sizeof x);

        x = this->geo.L6;
        out.write((char*)&x, sizeof x);
        x = this->geo.L7;
        out.write((char*)&x, sizeof x);

        x = this->geo.dd1;
        out.write((char*)&x, sizeof x);
        x = this->geo.dd2;
        out.write((char*)&x, sizeof x);
        x = this->geo.dd3;
        out.write((char*)&x, sizeof x);
        x = this->geo.dd4;
        out.write((char*)&x, sizeof x);
        x = this->geo.dd5;
        out.write((char*)&x, sizeof x);
        x = this->geo.dd6;
        out.write((char*)&x, sizeof x);
        x = this->geo.dd7;
        out.write((char*)&x, sizeof x);
        x = this->geo.dd8;
        out.write((char*)&x, sizeof x);
    }

    // Записываем количество узлов и все узлы
    n = this->All_Yzel.size();
    out.write((char*)&n, sizeof n);
    for (auto& i : this->All_Yzel)
    {
        x = i->coord[0][0];
        y = i->coord[0][1];
        out.write((char*)&x, sizeof x);
        out.write((char*)&y, sizeof y);
    }

    // Теперь записываем лучи (для движения сетки в 3Д программе) 
    cout << " this->All_name_luch.size() = " << this->All_name_luch.size() << endl;
    for (auto& LL : this->All_name_luch)
    {
        cout << "SSSSS  " << endl;
        auto L = *LL;
        n = L.size();
        cout << "n = " << n << endl;
        out.write((char*)&n, sizeof n);
        for (auto& i : L)
        {
            n = i->Yzels.size();
            out.write((char*)&n, sizeof n);
            for (auto& j : i->Yzels)
            {
                n = j->number;
                out.write((char*)&n, sizeof n);
            }

            n = i->Yzels_opor.size();
            out.write((char*)&n, sizeof n);
            cout << "n =   " << n <<  endl;
            for (auto& j : i->Yzels_opor)
            {
                n = j->number;
                out.write((char*)&n, sizeof n);
            }
        }
    }

    // Записываем количество ячеек и все ячейки (номера узлов, из которых состоит ячейка)
    n = this->All_Cell.size();
    out.write((char*)&n, sizeof n);
    for (auto& i : this->All_Cell)
    {
        n = i->Yzels.size();
        out.write((char*)&n, sizeof n);
        for (auto& j : i->Yzels)
        {
            n = j->number;
            out.write((char*)&n, sizeof n);
        }
    }

    // Теперь записываем всех соседей ячейки
    for (auto& i : this->All_Cell)
    {
        n = i->Soseds.size();
        out.write((char*)&n, sizeof n);
        for (auto& j : i->Soseds)
        {
            n = j->number;
            out.write((char*)&n, sizeof n);
        }
    }

    // нужны ли мне сорта ячеек? Кажется что это интересное знание, но, 
    // если например подразбивать сетку, то эти сорта становятся бесполезными


    int kkk = 1000;
    // Для будующих добавок какой-либо информации 
    // Если добавки нет, пишется ноль, если есть, пишется 1 и сама добавка
    // Это позволяет впоследствии считывать файл, даже если чего-то нет в середине вывода
    // (например, где-то будут пикапы, где-то нет)

    // Если есть добавка пишем так:
    // {---------------------------------
    //kkk--;
    //n = 1;
    //out.write((char*)&n, sizeof n);
    // Тело добавки
    // }---------------------------------

    n = 0;
    for (int i = 0; i < kkk; i++)
    {
        out.write((char*)&n, sizeof n);
    }


    out.close();
}

void Setka::Testing(void)
{
    // Тест 1
    // Для каждой ячейки проверим, что её набор точек образуют четырёхугольник в правильном порядке
    // Этот тест лучше смотреть руками при печати ячеек сетки, лень было тут реализовывать
}
