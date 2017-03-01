#ifndef SHAPETRACKER_H
#define SHAPETRACKER_H

#include <iostream>
#include <utility>
#include <vector>
#include <stdio.h>

#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "LS.h"
#include "CmdLine.h"

#include "construct_graph.h"
#include "mat_main.h"
#include "construct_graph_dt.h"
#include "BDSP_optimization.h"

using namespace std;
using namespace cv;

/// Function prototype for DetectEdgesByED exported by EDLinesLib.a
LS *DetectLinesByED(unsigned char *srcImg, int width, int height, int *pNoLines);

// Function for contructing graph on line segments extracted by EDLines
//vector<double*>* construct_graph_dt(vector<double*> Ltable, double WIDTH, double HEIGHT, double LAMBDA, double ALPHA, vector<Point> ptPolygon,vector<double*> &usedLines);

// Function for contour grouping by RC algorithm
//int RC(vector<double*> *graph, int ncount, int ecount, int NumIteration, string fname);

// Function for contour grouping by our algorithm: bidirectional shortest path searching
//int BDSP(vector<double*> *graph,vector<double*> lines_vc, vector<Point> &ptPolygon, int WIDTH, int HEIGHT);

int lines_detector(Mat grayScaleFrame, vector<double*> &lines_vc);

vector<int>* shp_tracker(Mat grayScaleframe, vector<double*> lines_vc, vector<Point> &shp_points, int gp_flag);

#endif // SHAPETRACKER_H
