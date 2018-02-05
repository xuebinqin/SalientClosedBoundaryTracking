#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <opencv2/opencv.hpp>

#include "Timer.h"
#include "LS.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <map>
#include "CmdLine.h"

#include <ctime>

#include <sstream>

#include "construct_graph.h"
#include "mat_main.h"
#include "construct_graph_dt.h"
#include "BDSP_optimization.h"
#include "shapeTracker.h"
#include "CommonFunctions.h"


using namespace cv;
using namespace std;

int video_proc = 0;//video_proc: 0 read from webcam, 1 read from video file
int gp_flag = 1;//0: RRC, 1: BDP

int k, pk=1;
int pki=1;
bool tracking_flg = false;

vector<Point> points; //contour of pre-shape
Point pt_lbd, pt_mv;
int mbd = 0;//flag of right button click

double ave_fps = 0;//average fps

void CallBackFunc(int event, int x, int y, int flags, void* userdata)
{
    if  ( event == EVENT_LBUTTONDOWN )
    {
        //Point* pt = (Point*) userdata;
        if(!tracking_flg){
            if(mbd == 0){
                pt_lbd.x = x;
                pt_lbd.y = y;
                points.push_back(pt_lbd);
            }
        }
        //cout << "Left button of the mouse is clicked - position (" << x << ", " << y << ")" << endl;
    }
    else if  ( event == EVENT_RBUTTONDOWN )
    {

        //cout << "Right button of the mouse is clicked - position (" << x << ", " << y << ")" << endl;
    }
    else if  ( event == EVENT_MBUTTONDOWN )
    {
        mbd = 1;//stop polygon selection
        //cout << "Middle button of the mouse is clicked - position (" << x << ", " << y << ")" << endl;
    }
    else if ( event == EVENT_MOUSEMOVE )
    {

        if(!tracking_flg){
            pt_mv.x = x;
            pt_mv.y = y;
        }

    }

    return;
}

int main(int argc, char* argv[]){

    if(argc==3||argc==4){
        video_proc = atoi(argv[1]);
        gp_flag = atoi(argv[2]);
    }
    else{
        cout<<"=====================WARNING======================"<<endl;
        cout<<"Please specify the parameters correctly as follows:"<<endl;
        cout<<"\"./SalientBoundaryTracking <video_proc> <gp_flag> <video_path>\""<<endl;
        cout<<"such as \"./SalientBoundaryTracking 1 1 ./TEST_VIDEO.avi\"."<<endl;
        cout<<"video_proc: "<<endl;
        cout<<"            \"0\": read file from your webcam"<<endl;
        cout<<"            \"1\": read video from a .avi file"<<endl;
        cout<<"gp_flag: "<<endl;
        cout<<"             \"0\": use RRC adapted method to track"<<endl;
        cout<<"             \"1\": use our BDSP to track"<<endl;
        cout<<"video_path: specify the video sequence need to be loaded"<<endl;
        cout<<"            if video_proc == 0, it will be neglected."<<endl;
        cout<<"-------------------------------------------------"<<endl;
        cout<<"(1). Click mouse left button to initialize a"<<endl;
        cout<<"     polygon in \"Display\" window."<<endl;
        cout<<"     When read video from a webcam,"<<endl;
        cout<<"     you have to press key \"i\" before"<<endl;
        cout<<"     initialization."<<endl;
        cout<<"(2). Click mouse middle button to finish initialization."<<endl;
        cout<<"(3). Press key \"space\" to run tracking continuously or"<<endl;
        cout<<"     press key \"s\" to run tracking frame by frame."<<endl;
        cout<<"     Key \"r\" can recover the"<<endl;
        cout<<"     continuously tracking model)."<<endl;
        cout<<"(4). Tracking output are two or three files: "<<endl;
        cout<<"     \"TRACKING_RESULTS.txt\""<<endl;
        cout<<"     \"RECORDED_RESULTS.avi\""<<endl;
        cout<<"     \"RECOREDED_VIDEO.avi\"(if video_proc==0)."<<endl;
        cout<<"-------------------------------------------------"<<endl;
        return 0;
    }

    clock_t start_line, start_tracking;
    double total_time_line = 0;
    double total_time_track = 0;
    double total_time = 0;

    double Graph_edges = 0;
    double Graph_nodes = 0;

    VideoCapture cap;

    if(video_proc==0){
        cap = VideoCapture(0);
    }
    else if(video_proc==1){
        const char* video_path = argv[3];
        cap = VideoCapture(video_path);
    }

    if(!cap.isOpened())  // check if we succeeded
        return -1;


    //Write Video_frame
    int frame_width=   cap.get(CV_CAP_PROP_FRAME_WIDTH);
    int frame_height=   cap.get(CV_CAP_PROP_FRAME_HEIGHT);

    VideoWriter video;
    if(video_proc==0){
        video = VideoWriter("RECORDED_VIDEO.avi",CV_FOURCC('M','J','P','G'),30,Size(frame_width,frame_height),true);
    }

    //Write Video_results

    VideoWriter videor;
    videor = VideoWriter("RECORDED_RESULTS.avi",CV_FOURCC('M','J','P','G'),30,Size(frame_width,frame_height),true);

    Mat frame;
    Mat grayscaleFrame;
    int width, height;

    int frame_whole_id=0;
    int frame_id=0;

    remove("TRACKING_RESULTS.txt");//remove prevoius results

    for(;;){

        //cap >> frame; // get a new frame from camera
        if(video_proc){
            if(frame_whole_id==0||tracking_flg){
                if(!cap.read(frame)){
                    break;
                }
            }
        }else{
            if(!cap.read(frame)){
                break;
            }
        }

        //- DETECT LINE SEGMENTS USING EDLINES-//
        cvtColor(frame,grayscaleFrame, CV_BGR2GRAY);//convert rgb image to gray image

        width = grayscaleFrame.size().width;
        height = grayscaleFrame.size().height;

        double line_time = 0;
        start_line = clock();
        //EDLines detector
        vector<double*> lines_vc;
        lines_detector(grayscaleFrame, lines_vc);

        line_time = (clock()-start_line)/(double)CLOCKS_PER_SEC;

//=================initialization==================
        if(!tracking_flg){

            Mat frame_ini;
            frame_ini = frame.clone();
            draw_lines(frame_ini, lines_vc,Scalar(0,0,255),1);//draw detected line segments


            if(points.size()>1){
                for(int i = 0;i<points.size()-1; i++){
                    line(frame_ini, points[i], points[i+1], Scalar(0,255,255), 2, CV_AA);
                }
            }

            if(mbd == 0){
                if(points.size()>0){
                    line(frame_ini, pt_mv, points[points.size()-1], Scalar(0,255,255), 2, CV_AA);
                    line(frame_ini, pt_mv, points[0], Scalar(0,255,255), 2, CV_AA);
                }
            }
            else
            {
                line(frame_ini, points[points.size()-1], points[0], Scalar(0,255,255), 2, CV_AA);
            }


            namedWindow("Display",1);
            imshow("Display",frame_ini);
            k = cvWaitKey(pki) & 255;

            //initialize the tracker
            if(k == 105||video_proc){//"i", fix the frame and initialize trackers
                //pki = 0;
                //set the callback function for any mouse event
                setMouseCallback("Display", CallBackFunc, &points);
            }

            //if(pki == 0 ){//space key, initialize tracker
            if(k == 32){//space key
                mbd = 0;
                tracking_flg = true;

                k = 1;
                pk = 1;

            }
            else if(k == 115){// s
                tracking_flg = true;
                k = 1;
                pk = 0;

            }
            //}
        }

        if(tracking_flg){
//================================start tracking======================================

            //start_frame = clock();
            cout<<"--- fame_id: "<<frame_id<<"---"<<endl;

            if(video_proc==0){//record video when read video stream from webcam
                video.write(frame);
            }

            draw_lines(frame, lines_vc,  Scalar(0,0,255), 1);//draw detected line segments
            //draw_contour(frame, points, Scalar(0,255,255), 2);//draw tracked contour of the last frame

            vector<int>* graphSize;

            double tracking_time = 0;
            start_tracking = clock();

            graphSize = shp_tracker(grayscaleFrame, lines_vc, points, gp_flag);//tracking "points" will be changed to contours of the current frame
            //cout<<"graph size, edge: "<< (*graphSize).at(0)<<" node: "<<(*graphSize).at(1)<<endl;

            tracking_time = (clock()-start_tracking)/(double)CLOCKS_PER_SEC;
            //cout<<"tracking time: "<<tracking_time<<" s"<<endl;

            //double frame_time;
            //frame_time = (clock()-start_tracking)/(double)CLOCKS_PER_SEC;
            //cout<<"frame time: "<<frame_time<<" s"<<endl;
            //cout<<"FPS: "<<1.0/frame_time<<endl;

            total_time_line += line_time;
            total_time_track  += tracking_time;


            total_time += line_time;
            total_time += tracking_time;

            Graph_edges += (*graphSize).at(0);
            Graph_nodes += (*graphSize).at(1);


            //cout<<"lines: "<<lines_vc.size();

            draw_contour(frame, points, Scalar(0,255,0), 1);//draw tracked contour of the current frame
            namedWindow("Display", WINDOW_AUTOSIZE);
            imshow("Display",frame);

//-------------save resulting edge map into a folder------------------
//            Mat rlt = Mat::ones(height, width, CV_8UC1);
//            for(int i = 0; i< points.size()-1; i++){
//                line(rlt, points[i], points[i+1],Scalar(255), 1, 8);
//            }
//            line(rlt, points[points.size()-1], points[0],Scalar(255), 1, 8);

//            char tmp_name[20];
//            sprintf(tmp_name,"%04d.png",frame_id+1);
//            imwrite("../../ShpTkr_data/RRCTracker_png/box_359/"+string(tmp_name),rlt);

/*
            char frmNameOverlap[20];
            sprintf(frmNameOverlap,"../SalientClosedBoundaryTracking_Results/BookStand/%04d.png",frame_id);
            string frmNameOPStr;
            frmNameOPStr = frmNameOverlap;
            imwrite(frmNameOPStr,frame);
*/
            save_track_results(points,  "TRACKING_RESULTS.txt");

            videor.write(frame);
            frame_id++;
            k = cvWaitKey(pk) & 255;

            if(k == 115){// "s"
                pk = 0;
            }
            else if(k == 114){// "r"
                pk = 1;
            }
            //break;
        }//if(!tracking_flag){}else{

        frame_whole_id++;

        if(k == 27){
            break;
            //return 0;
        }
        else if(k == 32){

        }

    }//for(;;)
    cout<<"============Summation============"<<endl;
    cout<<"total_frame: "<< frame_id<<endl;

    //cout<<"total edges: "<<Graph_edges<<endl;
    cout<<"average edges: "<<double(Graph_edges)/double(frame_id)<<endl;

    //cout<<"total nodes: "<<Graph_nodes<<endl;
    cout<<"average nodes: "<<double(Graph_nodes)/double(frame_id)<<endl;

    //cout<<"total line detection time: "<<total_time_line<<endl;
    cout<<"average line detection time (s): "<<double(total_time_line)/double(frame_id)<<endl;

    //cout<<"total track time: "<<total_time_track<<endl;
    cout<<"average tracking time (s): "<<double(total_time_track)/double(frame_id)<<endl;

    //cout<<"total_time: "<< total_time<<endl;
    cout<<"average fps: "<< (double)frame_id/total_time <<endl;

    return 0;

} //end-main
