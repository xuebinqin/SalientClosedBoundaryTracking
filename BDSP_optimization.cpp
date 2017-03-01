#include "BDSP_optimization.h"

// Function for contour grouping by our algorithm: bidirectional shortest path searching
int BDSP(vector<double*> *graph,vector<double*> lines_vc, vector<Point> &ptPolygon, int width, int height){

    int nlines = lines_vc.size();
    int nedges = (*graph).size();
    double lastArea = contourArea(ptPolygon);
    //double lastLength = arcLength(ptPolygon, true);

    //---------similarity measure: LAST FRAME-------//
    //draw prior polygon and distance transform
/*    Mat lastPoly_img, last_tmp, last_bw, lastDis_img;
    const Point* pptl[1] = {&ptPolygon[0]};
    int nptl[] = {(int)ptPolygon.size()};
    //image saves polygon of last frame
    lastPoly_img = Mat::zeros(height, width, CV_8UC1);
    polylines(lastPoly_img, pptl, nptl, 1, 1, Scalar(1), 1, 8, 0);
    //image saves distance transform of last frame contour
    last_tmp = 255 - lastPoly_img;
    threshold(last_tmp, last_bw, 254, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
    distanceTransform(last_bw, lastDis_img, CV_DIST_L2, 3);
*/    //-------------------------------------------------//

    typedef boost::property<boost::edge_weight_t, float> EdgeWeightProperty;

    typedef boost::adjacency_list < boost::listS, boost::vecS, boost::undirectedS,
            boost::no_property, EdgeWeightProperty > Graph;

    typedef boost::graph_traits < Graph >::vertex_descriptor vertex_descriptor;
    typedef boost::graph_traits < Graph >::edge_descriptor edge_descriptor;
    typedef std::pair<int, int> Edge;

    //Create a graph
    Graph g;

    vector<Graph::vertex_descriptor> Vtx;
    //add end points of line segments as vertexes
    //each line has two vertexes
    for(int i = 0; i < nlines; i++){
        Vtx.push_back(boost::add_vertex(g));
        Vtx.push_back(boost::add_vertex(g));
    }

    for(int i = 0; i < nedges;){
        EdgeWeightProperty weighti((*graph)[i][3]);

        int a = int((*graph)[i][0])%4;
        int b = int((*graph)[i][1])%4;
        int ai;
        int bi;
        if(a == 0 || a == 1){
            ai = int((*graph)[i][0])/4*2;
        }
        else if(a == 2 || a == 3){
            ai = int((*graph)[i][0])/4*2 + 1;
        }

        if(b == 0 || b == 1){
            bi = int((*graph)[i][1])/4*2;
        }
        else if(b == 2 || b == 3){
            bi = int((*graph)[i][1])/4*2 + 1;
        }

        boost::add_edge(Vtx[ai], Vtx[bi], weighti, g);//n=2*i+1
        //}
        i = i + 2;
    }

    vector<int> cycle_id_a,cycle_id_b,cycle_id;//ID of vertex
    vector<Point> cycle_pts;//For compute area
    double cycle_gap = 0;
    double cycle_area = 0;
    //double cycle_length = 0;
    double cycle_ratio = 0;


    for(int i = 0; i < nlines;){

        //set the edge weight of the current line to infinite
        boost::remove_edge(Vtx[2*i], Vtx[2*i+1],g);
        EdgeWeightProperty weightii(10000);
        boost::add_edge(Vtx[2*i], Vtx[2*i+1], weightii, g);


        //Bidirectional shortest path
        // Create things for Dijkstra
        std::vector<vertex_descriptor> parents_a(boost::num_vertices(g)); // To store parents
        std::vector<int> distances_a(boost::num_vertices(g)); // To store distances
        // Compute shortest paths from v(2i) to all vertices, and store the output in parents and distances
        boost::dijkstra_shortest_paths(g, Vtx[2*i], boost::predecessor_map(&parents_a[0]).distance_map(&distances_a[0]));

        // Create things for Dijkstra
        std::vector<vertex_descriptor> parents_b(boost::num_vertices(g)); // To store parents
        std::vector<int> distances_b(boost::num_vertices(g)); // To store distances
        // Compute shortest paths from v(2*i+1) to all vertices, and store the output in parents and distances
        boost::dijkstra_shortest_paths(g, Vtx[2*i+1], boost::predecessor_map(&parents_b[0]).distance_map(&distances_b[0]));

        //set the edge weight of the current line back to zero
        boost::remove_edge(Vtx[2*i], Vtx[2*i+1],g);
        EdgeWeightProperty weightiii(0);
        boost::add_edge(Vtx[2*i], Vtx[2*i+1], weightiii, g);


        for(int j = 0; j < 2*nlines; ){

            // Output results
            //cycle_id_a.push_back(2*i);
            //cycle_id_b.push_back(2*i+1);

            int end = j;
            int tmpa = end;
            int tmpb = end;
            cycle_gap = distances_a[end] + distances_b[end];//store the distance of cycle


            //store the third points first
            //cycle_id_a.push_back(end);
            //cycle_id.push_back(end);
            if(int(end)%2==0){
                cycle_pts.push_back(Point2f(lines_vc[int(end)/2][0], lines_vc[int(end)/2][1]));
            }
            else if(int(end)%2==1){
                cycle_pts.push_back(Point2f(lines_vc[int(end)/2][2], lines_vc[int(end)/2][3]));
            }

            if(j!=2*i){
                while(tmpa!=parents_a[tmpa]/*2*i*/){

                    //cycle_id_a.push_back(parents_a[tmpa]);
                    //cycle_id.push_back(parents_a[tmpa]);
                    //store corresponding points of a half cycle
                    if(int(parents_a[tmpa])%2==0){
                        cycle_pts.push_back(Point2f(lines_vc[int(parents_a[tmpa])/2][0], lines_vc[int(parents_a[tmpa])/2][1]));
                    }
                    else if(int(parents_a[tmpa])%2==1){
                        cycle_pts.push_back(Point2f(lines_vc[int(parents_a[tmpa])/2][2], lines_vc[int(parents_a[tmpa])/2][3]));
                    }

                    tmpa = parents_a[tmpa];
                }
            }

            //trace the other half cycle vertex
            while(tmpb!=parents_b[tmpb]/*2*i+1*/){
                cycle_id_b.push_back(parents_b[tmpb]);
                tmpb = parents_b[tmpb];
            }
            //cycle_id_b.push_back(end);

            //store corresponding points of another half cycle inverse order
            for(int t = cycle_id_b.size()-1; t > -1; t--){
                int tmpt = cycle_id_b[t];
                //cycle_id.push_back(tmpt);
                if(int(tmpt)%2==0){
                    cycle_pts.push_back(Point2f(lines_vc[int(tmpt)/2][0], lines_vc[int(tmpt)/2][1]));
                }
                else if(int(tmpt)%2==1){
                    cycle_pts.push_back(Point2f(lines_vc[int(tmpt)/2][2], lines_vc[int(tmpt)/2][3]));
                }
            }

            //cout<<"size of cycle_pts: "<< cycle_pts.size()<<endl;
            if(cycle_pts.size()>0){
                cycle_area = contourArea(cycle_pts);
                //cycle_length = arcLength(cycle_pts,true);
            }

            //cout<<"area ratio: "<<cycle_area/contourArea(ptPolygon)<<endl;

            if(cycle_area/lastArea>0.90 && lastArea/cycle_area>0.90/* && lastLength/cycle_length>0.90 && cycle_length/lastLength>0.90*/){

                //---------similarity measure CURRENT FRAME-------//
                //draw prior polygon and distance transform
 /*               Mat currentPoly_img;// current_tmp, current_bw, currentDis_img;
                const Point* pptc[1] = {&cycle_pts[0]};
                int nptc[] = {(int)cycle_pts.size()};
                //image saves polygon of last frame
                currentPoly_img = Mat::zeros(height, width, CV_8UC1);
                polylines(currentPoly_img, pptc, nptc, 1, 1, Scalar(1), 1, 8, 0);

                Mat col;
                lastDis_img.copyTo(col,currentPoly_img);
                Scalar Dis_col = sum(col);
                double distance = Dis_col[0];
*/
/*                //image saves distance transform of last frame contour
                current_tmp = 255 - currentPoly_img;
                threshold(current_tmp, current_bw, 254, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
                distanceTransform(current_bw, currentDis_img, CV_DIST_L2, 3);

                Mat col, loc;
                currentDis_img.copyTo(col, lastPoly_img);
                lastDis_img.copyTo(loc, currentPoly_img);

                Scalar Dis_col = sum(col);
                Scalar Dis_loc = sum(loc);

                double distance = (double)Dis_col[0];
                if(Dis_loc[0] > Dis_col[0]){
                    distance = Dis_loc[0];
                }
*/


                if(cycle_ratio==0){//indicate the first cycle
                    //cout<<"-----debug 1-----"<<endl;

                    cycle_ratio = cycle_gap/cycle_area;

                    ptPolygon.swap(cycle_pts);

                }
                else{


                    //-------------------------------------------------//

                    if(cycle_ratio > cycle_gap/cycle_area){

                        cycle_ratio = cycle_gap/cycle_area;

                        ptPolygon.swap(cycle_pts);

                    }
                }

            }



            cycle_id_a.clear();
            cycle_id_b.clear();
            cycle_id.clear();//ID of vertex
            cycle_pts.clear();//For compute area

            vector<int>().swap(cycle_id_a);
            vector<int>().swap(cycle_id_b);
            vector<int>().swap(cycle_id);
            vector<Point>().swap(cycle_pts);
            j = j + 2;
        }

        i = i + 2;
    }




    return ptPolygon.size();

}
