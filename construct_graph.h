#ifndef CONSTRUCT_GRAPH_H
#define CONSTRUCT_GRAPH_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <opencv2/opencv.hpp>//
#include <opencv2/imgproc/imgproc.hpp>//
#include <opencv2/highgui/highgui.hpp>//
#include <map>//

#include "CommonFunctions.h"


using namespace std;
using namespace cv;//

int dround(double);
double distance(double x1, double y1, double x2, double y2);
int minimal(double, double, double, double);
double minimum(double, double, double, double);
int sign(double);
double area(double P1x, double P1y, double P2x, double P2y);
double triangle_area(double Ax, double Ay, double Bx, double By, double Cx, double Cy);
bool quad_inside(double p1x,double p1y,double p2x,double p2y,double p3x,double p3y,double p4x,double p4y,double px,double py);

double curvature(double P0x, double P0y, double P1x, double P1y, double P2x, double P2y, double P3x, double P3y);
double* triangle_incenter(double Ax,double Ay,double Bx,double By,double Cx,double Cy);
double* find_intersection(double P1x, double P1y, double P2x, double P2y, double P3x, double P3y, double P4x, double P4y);
double* find_corner(double P1x, double P1y, double P2x, double P2y, double Ix, double Iy, vector<double*> Ctable);
double angle(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y);

double gradient(double P1x, double P1y, double P2x, double P2y, vector<double*> G);

int order_points(vector<double*> usedLines, vector<Point> &points, vector<double*> contourEdges);

vector<double*>* construct_graph(vector<double*> allLtable, vector<double*> Ctable, vector<double*> grad_img, double WIDTH, double LAMBDA, bool homogeneity, bool use_corners,/* string datafile,*/ double ALPHA, vector<Point> ptPolygon,vector<double*> &Ltable);

#endif // CONSTRUCT_GRAPH_H
