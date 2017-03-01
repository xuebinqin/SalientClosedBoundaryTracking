#ifndef BDSP_OPTIMIZATION_H
#define BDSP_OPTIMIZATION_H
#include <boost/config.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/properties.hpp>

#include <boost/property_map/property_map.hpp>

#include <iostream>
#include <utility>
#include <vector>
#include <stdio.h>

#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;
using namespace cv;

// Function for contour grouping by our algorithm: bidirectional shortest path searching
int BDSP(vector<double*> *graph,vector<double*> lines_vc, vector<Point> &ptPolygon, int width, int height);
#endif // BDSP_OPTIMIZATION_H
