#include "construct_graph_dt.h"
/*#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <map>

using namespace cv;
using namespace std;

const double MAX_DISTANCE = 50;
const double PI = 3.141592653589793;

//int dround(double);
double distance_dt(double x1, double y1, double x2, double y2);
//int minimal(double, double, double, double);
//double minimum(double, double, double, double);
//int sign(double);
double area_dt(double P1x, double P1y, double P2x, double P2y);
//double triangle_area(double Ax, double Ay, double Bx, double By, double Cx, double Cy);
//bool quad_inside(double p1x,double p1y,double p2x,double p2y,double p3x,double p3y,double p4x,double p4y,double px,double py);

//double curvature(double P0x, double P0y, double P1x, double P1y, double P2x, double P2y, double P3x, double P3y);
//double* triangle_incenter(double Ax,double Ay,double Bx,double By,double Cx,double Cy);
//double* find_intersection(double P1x, double P1y, double P2x, double P2y, double P3x, double P3y, double P4x, double P4y);
//double* find_corner(double P1x, double P1y, double P2x, double P2y, double Ix, double Iy, vector<double*> Ctable);
//double angle(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y);

//double gradient(double P1x, double P1y, double P2x, double P2y, vector<double*> G);

//int intersectionTest(const double *line1, const vector<double*> &Ltable);

struct lessLine{
    bool operator()(const int* lhs, const int* rhs) const{ //first compare x then compare y
        return (lhs[0] == rhs[0]) ? ((lhs[2] == rhs[2]) ? ((lhs[1] == rhs[1]) ? (lhs[3] < rhs[3]) : (lhs[1] < rhs[1])) : (lhs[2] < rhs[2])) : (lhs[0] < rhs[0]);
    }
};

struct lessXY{
    bool operator()(const int* lhs, const int*rhs) const{
        return (lhs[0] == rhs[0]) ? (lhs[1] < rhs[1]) : (lhs[0] < rhs[0]);
    }
};
*/
//int main(int argc, char **argv) {
//vector<double*>* construct_graph_dt(vector<double*> allLtable, vector<double*> Ctable, vector<double*> grad_img, double WIDTH, double HEIGHT, double LAMBDA, bool homogeneity, bool use_corners, string datafile, double ALPHA, vector<Point> ptPolygon, vector<double*> &Ltable)
vector<double*>* construct_graph_dt(vector<double*> allLtable, double WIDTH, double HEIGHT, double LAMBDA, double ALPHA, vector<Point> ptPolygon, vector<double*> &Ltable) {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Mat test(HEIGHT,WIDTH,CV_8UC3,Scalar(0,0,0));//CREATE TESTING IMAGE TO SHOW LINE SEGMENT EDGES AND THE DELAUNAY TRIANGULATION EDGES

    //double *dataread1, *dataread2;
    //double *dataout;
    double *dtdataout;//
    double dtMAX_W1 = 0, dtMAX_W2 = 0;//
    vector<double*> *Gtable = new vector<double*>;
    //double MAX_W1 = 0, MAX_W2 = 0;
    ofstream output;
    int PSI = 0;
    //if (homogeneity) PSI = 1;

    //if (use_corners) output.open(datafile.c_str());

    //Delaunay Triangulation by OpenCV
    Rect rect(0,0,WIDTH,HEIGHT);
    //rect = boundingRect(pre_shape);
    //rectangle(test, rect, Scalar(0,125,125),1, 8, 0);


    //Create an instance of Subdiv2D
    Subdiv2D subdiv(rect);
    //Delaunay Triangulation edges
    vector<Vec4f> edgeList;

    //Create a vector of points
    //Vector<Point2f> points;

    // PART 1 - SOLID EDGE CONSTRUCTION //

    map<int*, int, lessLine> maplns;
    int* sline;

    // for inserting in Subdiv2D
    vector<Point2f> points;

    map<int*, int, lessXY> mappts;
    int* spoint;
    int* epoint;

    //vector<double*> Ltable;
    double* oneOrgLine;
    vector<double> aveDists;//average distance of two end points and a middle point to the polygon edge
    vector<double> totalDists;//total distance of two end points and a middle point to the polygon edge
    vector<double> angles;//angles between line segments on current frame and contour tracked on last frame

    double pre_area = contourArea(ptPolygon);


    for(int i = 0; i<(int)allLtable.size();i++){

        double dist1 = pointPolygonTest( ptPolygon, Point2f(allLtable[i][0],allLtable[i][1]), true);
        double dist2 = pointPolygonTest( ptPolygon, Point2f(allLtable[i][2],allLtable[i][3]), true);
        double dist3 = pointPolygonTest( ptPolygon, Point2f((allLtable[i][0]+allLtable[i][2])/2,(allLtable[i][1]+allLtable[i][3])/2),true);
        double aveDist = (abs(dist1)+abs(dist2)+abs(dist3))/3.0;
        double totalDist = (abs(dist1) > abs(dist2))?(abs(dist1) > abs(dist3)?abs(dist1):abs(dist3)):(abs(dist2) > abs(dist3)?abs(dist2):abs(dist3));
        //double aveDist = totalDist;


        if(aveDist<30){//average distance to previous contour less than 20
            aveDists.push_back(aveDist);
            totalDists.push_back(totalDist);

            oneOrgLine = new double[4];
            oneOrgLine[0] = allLtable[i][0];
            oneOrgLine[1] = allLtable[i][1];
            oneOrgLine[2] = allLtable[i][2];
            oneOrgLine[3] = allLtable[i][3];
            Ltable.push_back(oneOrgLine);
        }

    }


    for (int i=0; i<(int)Ltable.size(); i++) {
        // oder the line segments coordinates by x and y

        if(Ltable[i][0]>Ltable[i][2]){
            float temp1 = Ltable[i][0];
            float temp2 = Ltable[i][1];
            Ltable[i][0] = Ltable[i][2];
            Ltable[i][1] = Ltable[i][3];
            Ltable[i][2] = temp1;
            Ltable[i][3] = temp2;
        }
        else if(Ltable[i][0]==Ltable[i][2]){
            if(Ltable[i][1]>Ltable[i][3]){
                float temp1 = Ltable[i][0];
                float temp2 = Ltable[i][1];
                Ltable[i][0] = Ltable[i][2];
                Ltable[i][1] = Ltable[i][3];
                Ltable[i][2] = temp1;
                Ltable[i][3] = temp2;
            }

        }

        //map of index i and line segment edge
        sline = new int[4];
        sline[0] = cvRound(Ltable[i][0]);
        sline[1] = cvRound(Ltable[i][1]);
        sline[2] = cvRound(Ltable[i][2]);
        sline[3] = cvRound(Ltable[i][3]);
        maplns[sline] = i;

        //map of i and end point
        spoint = new int[2];
        spoint[0] = cvRound(Ltable[i][0]);
        spoint[1] = cvRound(Ltable[i][1]);
        mappts[spoint] = 2*i;
        epoint = new int[2];
        epoint[0] = cvRound(Ltable[i][2]);
        epoint[1] = cvRound(Ltable[i][3]);
        mappts[epoint] = 2*i+1;


        points.push_back(Point2f(float(Ltable[i][0]),float(Ltable[i][1])));
        points.push_back(Point2f(float(Ltable[i][2]),float(Ltable[i][3])));

        //one real line segment generate two edges on the graph
        dtdataout = new double[5];
        dtdataout[0] = 4*i;
        dtdataout[1] = 4*i + 2;
        dtdataout[2] = 0;
        dtdataout[3] = 0;
        dtdataout[4] = 1; // solid
        (*Gtable).push_back(dtdataout);
        //cout<<"id1, id2, dtw1: "<<dtdataout[0]<<" "<<dtdataout[1]<<" "<<dtdataout[2]<<" "<<dtdataout[3]<<endl;//testing

        dtdataout = new double[5];
        dtdataout[0] = 4*i + 1;
        dtdataout[1] = 4*i + 3;
        dtdataout[2] = 0;
        dtdataout[3] = 0;
        dtdataout[4] = 1; // solid
        (*Gtable).push_back(dtdataout);

        //cout<<"id1, id2, dtw1: "<<dtdataout[0]<<" "<<dtdataout[1]<<" "<<dtdataout[2]<<" "<<dtdataout[3]<<endl;//testing


        ///////////////////////////////////////////////////////////////////////////////////////////////////
        //TESTING DRAW LINE SEGMENTS
        Point p1,p2;
        p1.x = cvRound(Ltable[i][0]);
        p1.y = cvRound(Ltable[i][1]);
        p2.x = cvRound(Ltable[i][2]);
        p2.y = cvRound(Ltable[i][3]);
        line(test,p1,p2,Scalar(0,0,255),1,8,0);
        //cout<<"red: "<<p1.x<<" "<<p1.y<<" "<<p2.x<<" "<<p2.y<<endl;//testing

    }//end-for order by x or y

    subdiv.insert(points);

    //////////////////////////////////////

    // PART 2 - DASHED EDGE CONSTRUCTION USING DELAUNAY TRIANGULATION //

    // PART 2-1 REMOVE OVERLAPED EDGES OVER DETECTED LINE SEGMENTS //
    vector<double*> dashedLinesDT;
    double *line1;
    int *dline;

    //GET DELAUNAY TRIANGULATION EDGES
    subdiv.getEdgeList(edgeList);

    //cout<<"Edges generated by DT before overlap filtering: "<<edgeList.size()<<endl;

    //int num_edge = 0;
    vector<Point> pt(2);
    for(size_t i = 0; i < edgeList.size(); i++){

        Vec4f t = edgeList[i];
        pt[0] = Point(cvRound(t[0]),cvRound(t[1]));
        pt[1] = Point(cvRound(t[2]),cvRound(t[3]));

        // order the line segments coordinates by x and y
        line1 = new double[4];
        if(t[0]>t[2]){
            line1[0] = t[2];
            line1[1] = t[3];
            line1[2] = t[0];
            line1[3] = t[1];
        }
        else if(t[0]==t[2]){
            if(t[1]>t[3]){
                line1[0] = t[2];
                line1[1] = t[3];
                line1[2] = t[0];
                line1[3] = t[1];
            }
            else{
                line1[0] = t[0];
                line1[1] = t[1];
                line1[2] = t[2];
                line1[3] = t[3];
            }
        }
        else{
            line1[0] = t[0];
            line1[1] = t[1];
            line1[2] = t[2];
            line1[3] = t[3];
        }

        dline = new int[4];
        dline[0] = cvRound(line1[0]);
        dline[1] = cvRound(line1[1]);
        dline[2] = cvRound(line1[2]);
        dline[3] = cvRound(line1[3]);

        if(rect.contains(pt[0]) && rect.contains(pt[1]) && maplns.count(dline)==0 /*&& (distance(line1[0],line1[1],line1[2],line1[3]) < MAX_DISTANCE)*/){

            //cout<<"gap distance: "<<distance(line1[0],line1[1],line1[2],line1[3])<<endl;//testing
            //int flag_int = intersectionTest(line1,Ltable);
            //cout<<"flag_int: "<<flag_int<<endl;
            //if(flag_int==0){
                dashedLinesDT.push_back(line1);//save delaunay edges which are not overlapped with detected line segments
                ////////////////////////////////////////////////////////////////////////////////////////
                //TESTING DRAW LINE SEGMENTS
                Point p3,p4;
                p3.x = cvRound(line1[0]);
                p3.y = cvRound(line1[1]);
                p4.x = cvRound(line1[2]);
                p4.y = cvRound(line1[3]);
                line(test,p3,p4,Scalar(0,255,0),1,8,0);
                //cout<<"green: "<<p3.x<<" "<<p3.y<<" "<<p4.x<<" "<<p4.y<<endl;//testing
            //}

        }


    }
    //DELAUNAY TRIANGULATION EDGES WITHOUT OVERLAPS WITH DETECTED LINE SEGMENT EDGES
    //cout<<"Edges generated by DT after overlap filtering: "<<dashedLinesDT.size()<<endl;

    //////////////////////////////////////////////////////////////////////////////////////////
    //imshow("Gap filling by Delaunay Triangulation",test);
    //cvWaitKey(0);


    // PART 2-2 PUT DASHED EDGES TO THE GRAPH
    int* one_pt1;
    int* one_pt2;
    double pt1[2],ptb1[2],pt2[2],ptb2[2];
    double dtw1, dtarea1, dtarea2, dtcarea, dtgaparea;//, alpha, beta;
    double dtG1 = 0, dtG2 = 0, dtgapG = 0;
    int dtn1, dtn2, dtcount_dashed = 0;
    double *dtC, dtcurv;//, *I, curv; dtM1[2], dtM2[2];

    double dtw1_h1,dtw1_h2;//distance to the polygon edge of solid edge end points and middle point

    for(int i = 0; i < dashedLinesDT.size(); i++){

        //first end point of a dashed line
        one_pt1 = new int[2];
        one_pt1[0] = cvRound(dashedLinesDT[i][0]);
        one_pt1[1] = cvRound(dashedLinesDT[i][1]);

        //search the solid line segments which are connected by these dashed edges
        if(mappts.count(one_pt1)==1){
            if(mappts[one_pt1]%2 == 0){
                pt1[0] = points[mappts[one_pt1]].x;
                pt1[1] = points[mappts[one_pt1]].y;
                ptb1[0] = points[mappts[one_pt1]+1].x;
                ptb1[1] = points[mappts[one_pt1]+1].y;

                dtw1_h1 = totalDists[mappts[one_pt1]/2]/2.0;
            }
            else{
                pt1[0] = points[mappts[one_pt1]].x;
                pt1[1] = points[mappts[one_pt1]].y;
                ptb1[0] = points[mappts[one_pt1]-1].x;
                ptb1[1] = points[mappts[one_pt1]-1].y;

                dtw1_h1 = totalDists[mappts[one_pt1]/2]/2.0;
            }
        }
        else{
            continue;
        }

        dtn1 = 2*(mappts[one_pt1]/2);
        dtarea1 = area_dt(pt1[0],pt1[1],ptb1[0],ptb1[1]);//area of one of the solid line which is connected to current dashed line
        //pt1 and ptb1 are two end points of a solid line segment
        //dtarea1 is the area which is enclosed by itself and its projection on x axis (trapezoid area)

        if (pt1[0] > ptb1[0])
            dtn1++;

        //second end point of the dashed line
        one_pt2 = new int[2];
        one_pt2[0] = cvRound(dashedLinesDT[i][2]);
        one_pt2[1] = cvRound(dashedLinesDT[i][3]);
        if(mappts.count(one_pt2)==1){
            if(mappts[one_pt2]%2 == 0){
                pt2[0] = points[mappts[one_pt2]].x;
                pt2[1] = points[mappts[one_pt2]].y;
                ptb2[0] = points[mappts[one_pt2]+1].x;
                ptb2[1] = points[mappts[one_pt2]+1].y;

                dtw1_h2 = totalDists[mappts[one_pt2]/2]/2.0;
            }
            else{
                pt2[0] = points[mappts[one_pt2]].x;
                pt2[1] = points[mappts[one_pt2]].y;
                ptb2[0] = points[mappts[one_pt2]-1].x;
                ptb2[1] = points[mappts[one_pt2]-1].y;

                dtw1_h2 = totalDists[mappts[one_pt2]/2]/2.0;
            }
        }
        else{
            continue;
        }
        dtn2 = 2*(mappts[one_pt2]/2);
        dtarea2 = area_dt(pt2[0],pt2[1],ptb2[0],ptb2[1]);//area of the other solid line which is connected to current dashed line

        if (pt2[0] > ptb2[0])
            dtn2++;

        //for each dashed line generate two corresponding graph edges

        dtC = new double[2];
        dtC[0] = -1; dtC[1] = -1;

        double dist1_dsh = pointPolygonTest( ptPolygon, Point2f(pt1[0],pt1[1]), true);
        double dist2_dsh = pointPolygonTest( ptPolygon, Point2f(pt2[0],pt2[1]), true);
        double dist3_dsh = pointPolygonTest( ptPolygon, Point2f((pt1[0]+pt2[0])/2,(pt1[1]+pt2[1])/2),true);
        double aveDist_dsh = (abs(dist1_dsh) > abs(dist2_dsh))?(abs(dist1_dsh) > abs(dist3_dsh)?abs(dist1_dsh):abs(dist3_dsh)):(abs(dist2_dsh) > abs(dist3_dsh)?abs(dist2_dsh):abs(dist3_dsh));
               aveDist_dsh = aveDist_dsh + dtw1_h1 + dtw1_h2;
        //dashed edges distance and the half distance of its connected two solid edges

        int flag_constrain = 0;//distance based shape constraint
        int flag_area = 0;// area based shape constraint

        // calculate W1: gap length, distance to the shape of last frame, angles between edges and that of the last frame
        if (dtC[0] < 0) {
            //dtM1[0] = (pt1[0]+ptb1[0])/2;
            //dtM2[0] = (pt2[0]+ptb2[0])/2;
            //dtM1[1] = (pt1[1]+ptb1[1])/2;
            //dtM2[1] = (pt2[1]+ptb2[1])/2;

            dtcurv = 0;
            if(flag_constrain==0){
                dtw1 = pow(distance_dt(pt1[0],pt1[1],pt2[0],pt2[1]),ALPHA) + LAMBDA*dtcurv;//w1 represents the gap length (numerator), LAMBDA: 0, ALPHA: 1
            }
            else if(flag_constrain==1){
                dtw1 = pow(distance_dt(pt1[0],pt1[1],pt2[0],pt2[1]),ALPHA) + aveDist_dsh + LAMBDA*dtcurv;//w1 represents the gap length (numerator), LAMBDA: 0, ALPHA: 1
            }
            //dtw1 = pow(distance(pt1[0],pt1[1],pt2[0],pt2[1]),ALPHA)*aveDist_dsh + LAMBDA*dtcurv;//w1 represents the gap length (numerator), LAMBDA: 0, ALPHA: 1
        } //else {

         //   if(flag_constrain==0){
          //      dtw1 = pow((distance(pt1[0],pt1[1],dtC[0],dtC[1]) + distance(dtC[0],dtC[1],pt2[0],pt2[1])), ALPHA);
          //  }else if(flag_constrain==1){
          //      dtw1 = pow((distance(pt1[0],pt1[1],dtC[0],dtC[1]) + distance(dtC[0],dtC[1],pt2[0],pt2[1])), ALPHA) + aveDist_dsh;
          //  }
            ///dtw1 = pow((distance(pt1[0],pt1[1],dtC[0],dtC[1]) + distance(dtC[0],dtC[1],pt2[0],pt2[1])), ALPHA)*aveDist_dsh;
        //}

        if(flag_area == 1){
            dtw1 = pre_area*dtw1;
        }

        // calculate W2
        dtgaparea = area_dt(pt1[0],pt1[1],pt2[0],pt2[1]);//area of the current dashed edge

        //if (dtC[0] < 0)
            dtcarea = 0;
        //else
        //    dtcarea = triangle_area(pt1[0],pt1[1],dtC[0],dtC[1],pt2[0],pt2[1]);

        if (pt1[0] > pt2[0]) {
            dtgaparea = -dtgaparea;
            dtcarea = -dtcarea;
            dtgapG = -dtgapG;
        }


        // create dashed edge
        if (pt1[0] > ptb1[0]) {

            dtdataout = new double[5];
            dtdataout[0] = 2*dtn1;

            if (pt2[0] > ptb2[0]) {
                dtdataout[1] = 2*dtn2 + 1;
                dtdataout[2] = (dtarea1 - dtarea2)/2 + dtgaparea + dtcarea - PSI*((dtG1 - dtG2)/2 + dtgapG);
                dtdataout[3] = dtw1;
                dtdataout[4] = 0; // dashed
                (*Gtable).push_back(dtdataout);

                 //cout<<"id1, id2, dtw1: "<<dtdataout[0]<<" "<<dtdataout[1]<<" "<<dtdataout[2]<<" "<<dtdataout[3]<<endl;//testing


                dtdataout = new double[5];
                dtdataout[0] = 2*dtn1 + 1;
                dtdataout[1] = 2*dtn2;
                dtdataout[2] = (-dtarea1 + dtarea2)/2 - dtgaparea - dtcarea - PSI*((-dtG1 + dtG2)/2 - dtgapG);
                dtdataout[3] = dtw1;
                dtdataout[4] = 0; // dashed
                (*Gtable).push_back(dtdataout);

                 //cout<<"id1, id2, dtw1: "<<dtdataout[0]<<" "<<dtdataout[1]<<" "<<dtdataout[2]<<" "<<dtdataout[3]<<endl;//testing

                //if (dtC[0] > 0) {
                 //   output << 2*dtn1 << " " << 2*dtn2+1 << " " << dtC[0] << " " << dtC[1] << endl;
                 //   output << 2*dtn1+1 << " " << 2*dtn2 << " " << dtC[0] << " " << dtC[1] << endl;
                //}

            }
            else {
                dtdataout[1] = 2*dtn2;
                dtdataout[2] = (dtarea1 + dtarea2)/2 + dtgaparea + dtcarea - PSI*((dtG1 + dtG2)/2 + dtgapG);
                dtdataout[3] = dtw1;
                dtdataout[4] = 0; // dashed
                (*Gtable).push_back(dtdataout);

                 //cout<<"id1, id2, dtw1: "<<dtdataout[0]<<" "<<dtdataout[1]<<" "<<dtdataout[2]<<" "<<dtdataout[3]<<endl;//testing

                dtdataout = new double[5];
                dtdataout[0] = 2*dtn1 + 1;
                dtdataout[1] = 2*dtn2 + 1;
                dtdataout[2] = (-dtarea1 - dtarea2)/2 - dtgaparea - dtcarea - PSI*((-dtG1 - dtG2)/2 - dtgapG);
                dtdataout[3] = dtw1;
                dtdataout[4] = 0; // dashed
                (*Gtable).push_back(dtdataout);

                 //cout<<"id1, id2, dtw1: "<<dtdataout[0]<<" "<<dtdataout[1]<<" "<<dtdataout[2]<<" "<<dtdataout[3]<<endl;//testing

                //if (dtC[0] > 0) {
                //    output << 2*dtn1 << " " << 2*dtn2 << " " << dtC[0] << " " << dtC[1] << endl;
                //   output << 2*dtn1+1 << " " << 2*dtn2+1 << " " << dtC[0] << " " << dtC[1] << endl;
                //}

            }

        } else {

            dtdataout = new double[5];
            dtdataout[0] = 2*dtn1 + 1;

            if (pt2[0] > ptb2[0]) {
                dtdataout[1] = 2*dtn2 + 1;
                dtdataout[2] = (-dtarea1 - dtarea2)/2 + dtgaparea + dtcarea - PSI*((-dtG1 - dtG2)/2 + dtgapG);
                dtdataout[3] = dtw1;
                dtdataout[4] = 0; // dashed
                (*Gtable).push_back(dtdataout);

                 //cout<<"id1, id2, dtw1: "<<dtdataout[0]<<" "<<dtdataout[1]<<" "<<dtdataout[2]<<" "<<dtdataout[3]<<endl;//testing

                dtdataout = new double[5];
                dtdataout[0] = 2*dtn1;
                dtdataout[1] = 2*dtn2;
                dtdataout[2] = (dtarea1 + dtarea2)/2 - dtgaparea - dtcarea - PSI*((dtG1 + dtG2)/2 - dtgapG);
                dtdataout[3] = dtw1;
                dtdataout[4] = 0; // dashed
                (*Gtable).push_back(dtdataout);

                //cout<<"id1, id2, dtw1: "<<dtdataout[0]<<" "<<dtdataout[1]<<" "<<dtdataout[2]<<" "<<dtdataout[3]<<endl;//testing

                //if (dtC[0] > 0) {
                //    output << 2*dtn1+1 << " " << 2*dtn2+1 << " " << dtC[0] << " " << dtC[1] << endl;
                //    output << 2*dtn1 << " " << 2*dtn2 << " " << dtC[0] << " " << dtC[1] << endl;
                //}

            }
            else {
                dtdataout[1] = 2*dtn2;
                dtdataout[2] = (-dtarea1 + dtarea2)/2 + dtgaparea + dtcarea - PSI*((-dtG1 + dtG2)/2 + dtgapG);
                dtdataout[3] = dtw1;
                dtdataout[4] = 0; // dashed
                (*Gtable).push_back(dtdataout);

                //cout<<"id1, id2, dtw1: "<<dtdataout[0]<<" "<<dtdataout[1]<<" "<<dtdataout[2]<<" "<<dtdataout[3]<<endl;//testing

                dtdataout = new double[5];
                dtdataout[0] = 2*dtn1;
                dtdataout[1] = 2*dtn2 + 1;
                dtdataout[2] = (dtarea1 - dtarea2)/2 - dtgaparea - dtcarea - PSI*((dtG1 - dtG2)/2 - dtgapG);
                dtdataout[3] = dtw1;
                dtdataout[4] = 0; // dashed
                (*Gtable).push_back(dtdataout);

                //cout<<"id1, id2, dtw1: "<<dtdataout[0]<<" "<<dtdataout[1]<<" "<<dtdataout[2]<<" "<<dtdataout[3]<<endl;//testing

                //if (dtC[0] > 0) {
                  //  output << 2*dtn1+1 << " " << 2*dtn2 << " " << dtC[0] << " " << dtC[1] << endl;
                  //  output << 2*dtn1 << " " << 2*dtn2+1 << " " << dtC[0] << " " << dtC[1] << endl;
                //}

            }

        }

        delete[] dtC;

        dtcount_dashed += 2;

        if (dtw1 > dtMAX_W1) dtMAX_W1 = dtw1;
        if (fabs(dtdataout[2]) > dtMAX_W2) dtMAX_W2 = fabs(dtdataout[2]);

    }


    //-------------------above: dash edges generation---------------------------

    //if (use_corners) output.close();

    //cout << "dtDashed edges: " << dtcount_dashed << endl;

    dtdataout = new double[5];
    dtdataout[0] = dtMAX_W2;
    dtdataout[1] = dtMAX_W1;
    dtdataout[2] = 0;
    dtdataout[3] = 0;
    dtdataout[4] = 1;
    (*Gtable).push_back(dtdataout);

    // << "dtMAX_W1: " << dtMAX_W1 << ", dtMAX_W2:" << dtMAX_W2 << endl;

    return Gtable;
}


double distance_dt(double x1, double y1, double x2, double y2) {
    return sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
}

double area_dt(double P1x, double P1y, double P2x, double P2y) {

    return ( fabs(P1x-P2x)*( (P1y>P2y ? P1y : P2y) - fabs(P1y-P2y)/2) );

}

/*
int minimal(double d1, double d2, double d3, double d4) {
    if (d1 <= d2 && d1 <= d3 && d1 <= d4){return 1;}
    if (d2 <= d1 && d2 <= d3 && d2 <= d4){return 2;}
    if (d3 <= d1 && d3 <= d2 && d3 <= d4){return 3;}
    if (d4 <= d1 && d4 <= d2 && d4 <= d3){return 4;}
    return 0;
}

double minimum(double d1, double d2, double d3, double d4) {
    if (d1 <= d2 && d1 <= d3 && d1 <= d4){return d1;}
    if (d2 <= d1 && d2 <= d3 && d2 <= d4){return d2;}
    if (d3 <= d1 && d3 <= d2 && d3 <= d4){return d3;}
    if (d4 <= d1 && d4 <= d2 && d4 <= d3){return d4;}
    return 0;
}

int sign(double d) {
    if (d == 0) return 0;
    return ((d>0)?1:-1);
}


int dround(double d) {
    return ((int)floor(d + 0.5));
}

double triangle_area(double Ax, double Ay, double Bx, double By, double Cx, double Cy) {

    return ( fabs(-Bx*Ay + Cx*Ay + Ax*By - Cx*By - Ax*Cy + Bx*Cy)/2 );

}



bool quad_inside(double p1x,double p1y,double p2x,double p2y,double p3x,double p3y,double p4x,double p4y,double px,double py) {

    int s1, s2, s3;

    // first, check if p inside triangle p1p2p3

    s1 = sign(-p2x*p1y + px*p1y + p1x*p2y - px*p2y - p1x*py + p2x*py);
    s2 = sign(-p3x*p2y + px*p2y + p2x*p3y - px*p3y - p2x*py + p3x*py);
    s3 = sign(-p1x*p3y + px*p3y + p3x*p1y - px*p1y - p3x*py + p1x*py);

    if ((s1 == s2) && (s1 == s3)) { // point is inside
        return true;
    } else { // check second triangle, p3p4p1

        s1 = sign(-p3x*p1y + px*p1y + p1x*p3y - px*p3y - p1x*py + p3x*py);
        s2 = sign(-p4x*p3y + px*p3y + p3x*p4y - px*p4y - p3x*py + p4x*py);
        s3 = sign(-p1x*p4y + px*p4y + p4x*p1y - px*p1y - p4x*py + p1x*py);

        if ((s1 == s2) && (s1 == s3))
            return true;
        else
            return false;

    }

}

double gradient(double P1x, double P1y, double P2x, double P2y, vector<double*> G) {

    double sum = 0, *line, m, b;
    int minx, maxx, maxy;

    if (P1x < P2x) {
        minx = dround(P1x);
        maxx = dround(P2x);
    } else {
        minx = dround(P2x);
        maxx = dround(P1x);
    }

    if (minx == maxx) return sum;

    m = (P1y-P2y)/(P1x-P2x);
    b = P1y - m*P1x;

    for (int i=minx; i<=maxx; i++) {

        maxy = dround(m*i + b);

        for (int j=0; j <= maxy; j++) {
            line = G[j];
            sum += line[i];
        }

    }

    return sum;

}

double* find_intersection(double P1x, double P1y, double P2x, double P2y, double P3x, double P3y, double P4x, double P4y) {

    double *I = new double[2];
    double m1,m2,b1,b2;

    if (P1x == P2x)
        m1 = 0;
    else {
        m1 = (P1y - P2y)/(P1x - P2x);
        b1 = P1y - m1*P1x;
    }
    if (P3x == P4x)
        m2 = 0;
    else {
        m2 = (P3y - P4y)/(P3x - P4x);
        b2 = P3y - m2*P3x;
    }

    if (P1x == P2x) {
        if (P3x == P4x) { //parallel, no intersection
            I[0] = -1;
            I[1] = -1;
        } else {
            I[0] = P1x;
            I[1] = m2*I[0] + b2;
        }
    } else {
        if (P3x == P4x) {
            I[0] = P3x;
            I[1] = m1*I[0] + b1;
        } else {
            if (m1 == m2) { //parallel, no intersection
                I[0] = -1;
                I[1] = -1;double area(double P1x, double P1y, double P2x, double P2y) {

    return ( fabs(P1x-P2x)*( (P1y>P2y ? P1y : P2y) - fabs(P1y-P2y)/2) );

}
            } else {
                I[0] = (b2-b1)/(m1-m2);
                I[1] = m1*I[0] + b1;
            }
        }
    }
    for (int i=0; i<(int)Ltable.size(); i++) {
        // oder the line segments coordinates by x and y

        if(Ltable[i][0]>Ltable[i][2]){
            float temp1 = Ltable[i][0];
            float temp2 = Ltable[i][1];
            Ltable[i][0] = Ltable[i][2];
            Ltable[i][1] = Ltable[i][3];
            Ltable[i][2] = temp1;
            Ltable[i][3] = temp2;
        }
        else if(Ltable[i][0]==Ltable[i][2]){
            if(Ltable[i][1]>Ltable[i][3]){
                float temp1 = Ltable[i][0];
                float temp2 = Ltable[i][1];
                Ltable[i][0] = Ltable[i][2];
                Ltable[i][1] = Ltable[i][3];
                Ltable[i][2] = temp1;
                Ltable[i][3] = temp2;
            }

        }

        //map of index i and line segment edge
        sline = new int[4];
        sline[0] = cvRound(Ltable[i][0]);
        sline[1] = cvRound(Ltable[i][1]);
        sline[2] = cvRound(Ltable[i][2]);
        sline[3] = cvRound(Ltable[i][3]);
        maplns[sline] = i;

        //map of i and end point
        spoint = new int[2];
    // discard I if P1,P2 and I are aligned
    if (((I[0] - P1x) < 0.01) && ((I[1] - P1y) < 0.01)) {
        I[0] = -1; I[1] = -1;
    } else if (((I[0] - P3x) < 0.01) && ((I[1] - P3y) < 0.01)) {
        I[0] = -1; I[1] = -1;
    }

    return I;
}

double* find_corner(double P1x, double P1y, double P2x, double P2y, double Ix, double Iy, vector<double*> Ctable) {

    double *C = new double[2], *p;
    double c1[2],c2[2],c3[2],c4[2];
    double m,b1,b2, d = 0;

    C[0] = -1;
    C[1] = -1;

    if (Ix == -1) {
        return C;
    }

    if (distance(P1x,P1y,Ix,Iy) > distance(P2x,P2y,Ix,Iy)) {
        c1[0] = P1x;
        c1[1] = P1y;
    } else {
        c1[0] = P2x;
        c1[1] = P2y;
    }

    double *inctr = triangle_incenter(P1x,P1y,Ix,Iy,P2x,P2y);

    if (Ix == inctr[0]) { //angle bisector is vertical
        c2[0] = 2*Ix - c1[0];
        c2[1] = c1[1];
    } else {
        m = (Iy-inctr[1])/(Ix-inctr[0]);
        if (m == 0) {
            c2[0] = c1[0];
            c2[1] = 2*Iy - c1[1];
        } else {
            b1 = c1[1] + c1[0]/m;
            b2 = Iy - m*Ix;
            c2[0] = (b2-b1)/(-1/m - m);
            c2[1] = m*c2[0] + b2;
        }
    }

    c3[0] = 2*Ix - c1[0];
    c3[1] = 2*Iy - c1[1];
    c4[0] = 2*Ix - c2[0];
    c4[1] = 2*Iy - c2[1];
*/
    /*if (P1x == 84 && P2x == 84) {
cout << "quad" << endl;
cout << c1[0] << " " << c1[1] << endl;
cout << c2[0] << " " << c2[1] << endl;
cout << c3[0] << " " << c3[1] << endl;
cout << c4[0] << " " << c4[1] << endl;
}*/

    //Search area is {c1,c2,c3,c4}
/*
    for (int i=0; i<(int)Ctable.size(); i++) {

        p = Ctable[i];

        if (p[2] > d) {
            if (quad_inside(c1[0],c1[1],c2[0],c2[1],c3[0],c3[1],c4[0],c4[1],p[0],p[1])) {
                C[0] = p[0];
                C[1] = p[1];
                d = p[2];
            }
        }
    }

    return C;

}


double* triangle_incenter(double Ax,double Ay,double Bx,double By,double Cx,double Cy) {

    double a,b,c, perimeter;
    double *center = new double[2];

    a = distance(Ax,Ay,Bx,By);
    b = distance(Bx,By,Cx,Cy);
    c = distance(Cx,Cy,Ax,Ay);

    perimeter = a+b+c;

    if ( perimeter == 0) {
        center[0] = Ax;
        center[1] = Ay;
    } else {
        center[0] = (b*Ax + c*Bx + a*Cx)/perimeter;
        center[1] = (b*Ay + c*By + a*Cy)/perimeter;
    }

    return center;
}


double curvature(double P0x, double P0y, double P1x, double P1y, double P2x, double P2y, double P3x, double P3y) {

    double total = 0, delta, x, xx, y, yy, K;
    double interm[4][2];
    int n = 100;
    double length, u, x1, y1, x2, y2;

    delta = 1.0 / n;

    interm[0][0] = -(P0x) + (3 * P1x) - (3 * P2x) + P3x;
    interm[0][1] = -(P0y) + (3 * P1y) - (3 * P2y) + P3y;
    interm[1][0] = (3 * P0x) - (6 * P1x) + (3 * P2x);
    interm[1][1] = (3 * P0y) - (6 * P1y) + (3 * P2y);
    interm[2][0] = -(3 * P0x) + (3 * P1x);
    interm[2][1] = -(3 * P0y) + (3 * P1y);
    interm[3][0] = P0x;
    interm[3][1] = P0y;

    for (int i = 0; i < n; i++)
    {
        x = 3 * (delta * delta * i * i) * interm[0][0] + 2 * delta * i * interm[1][0] + interm[2][0];
        xx = 6 * delta * i * interm[0][0] + 2 * interm[1][0];
        y = 3 * (delta * delta * i * i) * interm[0][1] + 2 * delta * i * interm[1][1] + interm[2][1];
        yy = 6 * delta * i * interm[0][1] + 2 * interm[1][1];

        if ( (x == 0) && (y == 0) )
            K = 0;
        else
        {
            K = ( fabs(x * yy - xx * y) / pow((fabs(x * x + y * y)), 1.5) );
            K = K * K;
        }

        u = delta*i;
        x1 = interm[0][0] * u * u * u + interm[1][0] * u * u + interm[2][0] * u + interm[3][0];
        y1 = interm[0][1] * u * u * u + interm[1][1] * u * u + interm[2][1] * u + interm[3][1];
        u = delta*(i+1);
        x2 = interm[0][0] * u * u * u + interm[1][0] * u * u + interm[2][0] * u + interm[3][0];
        y2 = interm[0][1] * u * u * u + interm[1][1] * u * u + interm[2][1] * u + interm[3][1];

        length = distance(x1,y1,x2,y2);

        total += K * length;
    }

    return total;

}

double angle(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y){
    double p1,p2, value;

    p1 = (p3x-p2x)*(p1x-p2x) + (p3y-p2y)*(p1y-p2y);
    p2 = (p3x-p2x)*(p1y-p2y) - (p3y-p2y)*(p1x-p2x);

    if (p1 == 0 && p2 == 0) return 0;

    value = atan2(p2,p1);

    return value;
}

int intersectionTest(const double *line1, const vector<double*> &Ltable){

intersectionTest
    //cout<<"line1: "<<line1[0]<<" "<<line1[1]<<" "<<line1[2]<<" "<<line1[3]<<" "<<endl;
    double p11[3] = {line1[0],line1[1],1};
    double p12[3] = {line1[2],line1[3],1};
    //the parameter of line1: a b c
    Mat P11(1,3,CV_64FC1,p11);
    Mat P12(1,3,CV_64FC1,p12);

    //cout<<"L1 points1: "<<P11.at<double>(0,0)<<" "<<P11.at<double>(0,1)<<" "<<P11.at<double>(0,2)<<endl;
    //cout<<"L1 points2: "<<P12.at<double>(0,0)<<" "<<P12.at<double>(0,1)<<" "<<P12.at<double>(0,2)<<endl;

    Mat L1 = P11.cross(P12);
    L1.at<double>(0,0) = L1.at<double>(0,0)/L1.at<double>(0,2);
    L1.at<double>(0,1) = L1.at<double>(0,1)/L1.at<double>(0,2);
    L1.at<double>(0,2) = 1;
    //cout<<"L1: "<<L1.at<double>(0,0)<<" "<<L1.at<double>(0,1)<<" "<<L1.at<double>(0,2)<<endl;

    int flag = 0;
    for(int i=0;i<Ltable.size();i++){
        double p21[3] = {Ltable[i][0],Ltable[i][1],1};
        double p22[3] = {Ltable[i][2],Ltable[i][3],1};

        Mat P21(1,3,CV_64FC1,p21);
        Mat P22(1,3,CV_64FC1,p22);

        //cout<<"L2 points1: "<<P21.at<double>(0,0)<<" "<<P21.at<double>(0,1)<<" "<<P21.at<double>(0,2)<<endl;
        //cout<<"L2 points2: "<<P22.at<double>(0,0)<<" "<<P22.at<double>(0,1)<<" "<<P22.at<double>(0,2)<<endl;

        Mat L2 = P21.cross(P22);
        L2.at<double>(0,0) = L2.at<double>(0,0)/L2.at<double>(0,2);
        L2.at<double>(0,1) = L2.at<double>(0,1)/L2.at<double>(0,2);
        L2.at<double>(0,2) = 1;

        //cout<<"L2: "<<L2.at<double>(0,0)<<" "<<L2.at<double>(0,1)<<" "<<L2.at<double>(0,2)<<endl;

        Mat Pt = L1.cross(L2);
        Pt.at<double>(0,0) = Pt.at<double>(0,0)/Pt.at<double>(0,2);
        Pt.at<double>(0,1) = Pt.at<double>(0,1)/Pt.at<double>(0,2);
        Pt.at<double>(0,2) = 1;

        double deltax1 = (Pt.at<double>(0,0)-P11.at<double>(0,0))*(Pt.at<double>(0,0)-P12.at<double>(0,0));
        double deltay1 = (Pt.at<double>(0,1)-P11.at<double>(0,1))*(Pt.at<double>(0,1)-P12.at<double>(0,1));

        double deltax2 = (Pt.at<double>(0,0)-P21.at<double>(0,0))*(Pt.at<double>(0,0)-P22.at<double>(0,0));
        double deltay2 = (Pt.at<double>(0,1)-P21.at<double>(0,1))*(Pt.at<double>(0,1)-P22.at<double>(0,1));

        //cout<<"deltax1: "<<deltax1<<endl;
        //cout<<"deltay1: "<<deltay1<<endl;

        //cout<<"deltax2: "<<deltax2<<endl;
        //cout<<"deltay2: "<<deltay2<<endl;

        //cout<<"intersection points: "<<Pt.at<double>(0,0)<<" "<<Pt.at<double>(0,1)<<" "<<Pt.at<double>(0,2)<<endl;

        //cvWaitKey(0);

        double epsilon = 1e-3;

        if(abs(Pt.at<double>(0,0)-P11.at<double>(0,0))>epsilon && abs(Pt.at<double>(0,0)-P12.at<double>(0,0))>epsilon &&
          abs(Pt.at<double>(0,1)-P11.at<double>(0,1))>epsilon && abs(Pt.at<double>(0,1)-P12.at<double>(0,1))>epsilon &&
          abs(Pt.at<double>(0,0)-P21.at<double>(0,0))>epsilon && abs(Pt.at<double>(0,0)-P22.at<double>(0,0))>epsilon &&
                abs(Pt.at<double>(0,1)-P21.at<double>(0,1))>epsilon && abs(Pt.at<double>(0,1)-P22.at<double>(0,1))>epsilon){
            if((deltax1<0||deltay1<0)&&(deltax2<0||deltay2<0)){
                flag = 1;
                cout<<"L1"<<endl;
                cout<<(Pt.at<double>(0,0)-P11.at<double>(0,0))<<endl;
                cout<<(Pt.at<double>(0,0)-P12.at<double>(0,0))<<endl;
                cout<<(Pt.at<double>(0,1)-P11.at<double>(0,1))<<endl;
                cout<<(Pt.at<double>(0,1)-P12.at<double>(0,1))<<endl;
                cout<<"L2"<<endl;
                cout<<(Pt.at<double>(0,0)-P21.at<double>(0,0))<<endl;
                cout<<(Pt.at<double>(0,0)-P22.at<double>(0,0))<<endl;
                cout<<(Pt.at<double>(0,1)-P21.at<double>(0,1))<<endl;
                cout<<(Pt.at<double>(0,1)-P22.at<double>(0,1))<<endl;
            }
        }

    }
    return flag;
}
*/
