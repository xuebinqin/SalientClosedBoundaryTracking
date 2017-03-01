#include "shapeTracker.h"
/*#include <iostream>
#include <utility>
#include <vector>
#include <stdio.h>

#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "LS.h"
#include "CmdLine.h"

using namespace std;
using namespace cv;

/// Function prototype for DetectEdgesByED exported by EDLinesLib.a
LS *DetectLinesByED(unsigned char *srcImg, int width, int height, int *pNoLines);

// Function for contructing graph on line segments extracted by EDLines
vector<double*>* construct_graph_dt(vector<double*> Ltable, double WIDTH, double HEIGHT, double LAMBDA, double ALPHA, vector<Point> ptPolygon,vector<double*> &usedLines);

// Function for contour grouping by RC algorithm
int RC(vector<double*> *graph, int ncount, int ecount, int NumIteration, string fname);

// Function for contour grouping by our algorithm: bidirectional shortest path searching
int BDSP(vector<double*> *graph,vector<double*> lines_vc, vector<Point> &ptPolygon, int WIDTH, int HEIGHT);
*/

int lines_detector(Mat grayScaleFrame, vector<double*> &lines_vc){

    vector<double*> tmplines;
    int width = grayScaleFrame.size().width;
    int height = grayScaleFrame.size().height;

    unsigned char *srcImg;
    srcImg = grayScaleFrame.data;

    int noLines;
    LS *lines = DetectLinesByED(srcImg, width, height, &noLines);

    double *one_line;

    if(noLines > 0){
        for (int i=0; i<noLines; i++){
            one_line = new double[4];
            one_line[0] = lines[i].sx;
            one_line[1] = lines[i].sy;
            one_line[2] = lines[i].ex;
            one_line[3] = lines[i].ey;
            tmplines.push_back(one_line);
        } //end-for
    }
    lines_vc.swap(tmplines);
    delete lines;
    return noLines;
}

//
vector<int>*  shp_tracker(Mat grayScaleframe, vector<double*> lines_vc, vector<Point> &shp_points, int gp_flag){

    int width = grayScaleframe.size().width;
    int height = grayScaleframe.size().height;

    //detected lines n*4
    //lines_detector(grayScaleframe, lines_vc);//lines detector
    vector<double*> *graph;//graph: node1 node2 area gap_length dash/solid xxx
    vector<double*> usedLines;//filtered lines m*4

    vector<double*> corners;
    vector<double*> gradient;

    if(gp_flag==0){
        graph = construct_graph(lines_vc, corners, gradient, width, 0.0, false, false, 1.0, shp_points, usedLines);
    }
    else{
        graph = construct_graph_dt(lines_vc, width, height, 0.0, 1.0, shp_points, usedLines);
    }


    double MAX_W2, MAX_W1;
    double *dataout;
    //vector<Point> plyg_points;//store the end points of grouped line segments subsequently as a polygon

    //get the last edge which is the initilization of iterations
    dataout = (*graph)[(*graph).size()-1];
    MAX_W2 = dataout[0];
    MAX_W1 = dataout[1];

    (*graph).pop_back();
    delete[] dataout;

    for (int i = 0; i < (int)(*graph).size(); i++) {
        dataout = (*graph)[i];
        dataout[2] = dataout[2]*600.0/MAX_W2;//normalize the area term
        dataout[3] = dataout[3]*600.0/MAX_W1;//normalize the gap length term
    }

    if(gp_flag==0){
        //RC(graph, (4*usedLines.size()), (*graph).size(), 1, fname);
        vector<double*> groupedEdges;
        RC(graph, (4*usedLines.size()), (*graph).size(), 1, groupedEdges);
        order_points(usedLines, shp_points, groupedEdges);
    }
    else{

        BDSP(graph, usedLines, shp_points, width, height);
    }
    //graph size returning
    vector<int> *graphSize = new vector<int>;
    (*graphSize).push_back((*graph).size());
    (*graphSize).push_back(2*usedLines.size());
    //cout<<"UsedLines number: "<<usedLines.size()<<endl;
    return graphSize;
}
