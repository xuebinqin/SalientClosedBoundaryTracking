#include "CommonFunctions.h"

int draw_lines(Mat frame, vector<double*> lines_vc, Scalar line_color, int line_width){
    if(lines_vc.size()>0){
        for(int i = 0; i < lines_vc.size(); i++){
            Point pt1, pt2;
            pt1.x = lines_vc[i][0];
            pt1.y = lines_vc[i][1];
            pt2.x = lines_vc[i][2];
            pt2.y = lines_vc[i][3];
            line(frame, pt1, pt2, line_color, line_width, CV_AA);
        }
    }
    else{
        cout<<"draw_lines failed: the size of lines_vc is zero!"<<endl;
    }

    return 1;
}

int draw_contour(Mat frame, vector<Point> points, Scalar point_color, int line_width){
    if(points.size()>0){
        for(int i = 0; i< points.size()-1; i++){
            line(frame, points[i], points[i+1],point_color, line_width, CV_AA);
        }
        line(frame, points[points.size()-1], points[0],point_color, line_width, CV_AA);
    }
    else{
        cout<<"draw_contour failed: the size of points is zero!"<<endl;
    }

    return 1;
}


int save_track_results(vector<Point> points, string fname){

    const char* filename = fname.c_str();

    FILE *fpt = fopen(filename, "a");
    fprintf(fpt, "%d ", points.size());
    for (int i=0; i<points.size(); i++){
        fprintf(fpt, "%6.2lf %6.2lf ", (float)points[i].x,(float)points[i].y);
    }
    fprintf(fpt,"\n");
    fclose(fpt);

    return 1;
}
