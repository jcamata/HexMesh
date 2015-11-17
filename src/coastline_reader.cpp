#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include  <cstdio>
using namespace std;

#include "coastline.h"

void CoastLineSetUp(CoastLine *cl)
{
    cl->npoints = 0;
    cl->pmin.x = cl->pmin.y = 1.0E+16;
    cl->pmax.x = cl->pmax.y = -1.0E+16;
}


void SetBBox(Point2D *pmin, Point2D *pmax, Point2D *p)
{
    pmin->x = (p->x < pmin->x)?p->x:pmin->x;
    pmin->y = (p->y < pmin->y)?p->y:pmin->y;
    pmax->x = (p->x > pmax->x)?p->x:pmax->x;
    pmax->y = (p->y > pmax->y)?p->y:pmax->y;
    
}

void CoastLineSetUpBBox(CoastLine *cl, Point2D* p)
{
    SetBBox(&cl->pmin, &cl->pmax, p);    
}

void Coastline_Reader(const char* file, SetOfCoastLine* scl)
{
    string x , y;
    fstream             in;
    in.open(file, ios::in);
    if(in.fail())
    {
        cerr << "ERROR: The file " << file<< " can not be openned!" << endl;
        exit(1);
    }
 scl->npart = 0;
    int n     = -1;
    while( in >> x >> y ) {

        //cout << x << "   " << y << endl; 

        if (x == "nan" && y == "nan" )  // new coastline
        {
            CoastLine line;
            CoastLineSetUp(&line);
            scl->coastlines.push_back(line);
            scl->npart++;
            n++;
            continue;
        }

        Point2D p;
        p.x = atof(x.c_str());
        p.y = atof(y.c_str());
        scl->coastlines[n].points.push_back(p);
        CoastLineSetUpBBox(&scl->coastlines[n], &p);
    }

    //cout << "It was found " << scl->npart <<" coastlines " <<endl;
    in.close();
 for(int i=0; i < scl->npart; i++)
    {
        scl->coastlines[i].npoints = scl->coastlines[i].points.size();
        //cout<< "Costline #"<<i<< " has " << scl->coastlines[i].npoints << " points" <<endl;
        //cout<< "  Bounding box: ("<<scl->coastlines[i].pmin.x <<" , "<<scl->coastlines[i].pmin.y << ") -- "
        //                     <<"("<<scl->coastlines[i].pmax.x <<" , "<<scl->coastlines[i].pmax.y << ")" << endl;
    }
}

