/* 
 * File:   coastline.h
 * Author: lucioc
 *
 * Created on November 16, 2015, 12:49 PM
 */

#ifndef COASTLINE_H
#define	COASTLINE_H

#include <stdio.h>
#include <stdlib.h>
#include <vector>

typedef struct {
    double x;
    double y;
} Point2D;

typedef struct {
    vector <Point2D> points;
    int             npoints;
    Point2D         pmin;
    Point2D         pmax;
} CoastLine;

typedef struct {
    vector<CoastLine>  coastlines;
    int                npart;
} SetOfCoastLine;

#endif	/* COASTLINE_H */

