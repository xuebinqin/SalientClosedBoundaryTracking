Slient Closed Boundary Tracking
============================
This code is implemeted for the papper `"Slaient Closed Boundary Tracking via Line Segments Perceptual Grouping, Xuebin Qin, Shida He, Camilo Perez Quintero, Abhineet Singh, Masood Dehghan and Martin Jagersand."` which is submited to IROS 2017.

Abstract
-----------------------------------------------
This paper presents a novel real-time method for tracking salient closed boundaries from video image sequences. This method operates on a set of straight line segments that are produced by line detection. The tracking scheme is coherently integrated into a perceptual grouping framework in which the visual tracking problem is tackled by identifying a subset of these line segments and connecting them sequentially to form a closed boundary with the largest saliency and a certain similarity to the previous one. Specifically, we define a new tracking criteria which combines a grouping cost and an area similarity constraint. This criteria makes the resulting boundary tracking more robust to local minima. To achieve real-time tracking performance, we use Delaunay Triangulation to build a graph model with the detected line segments and then reduce the tracking problem to finding the optimal cycle in this graph. This is solved by our newly proposed closed boundary candidates searching algorithm called "Bidirectional Shortest Path (BDSP)". The efficiency and robustness of the proposed method are tested on real video sequences as well as during a robot arm pouring experiment.

Used libraries
----------------------------------------------
This implementation are tested on `Opencv 2.4.9(3.1.0)`, `Boost 1.63.0` and `ubuntu 14.04 64 bit`.

The line detector used here is `EDLines(EDLines.a)` which is proposed in <br>
	`"C. Akinlar and C. Topal, EDLines: A real-time line segment detector with a false detection control, Pattern Recognition Letters, vol. 32, no. 13, pp. 1633-1642, 2011."` <br>
and can be downloaded from http://ceng.anadolu.edu.tr/CV/downloads/downloads.aspx. We include the lib in the root directory. It is worthynote that we are using the 64 bit ubuntu version of EDLines. We sugggest to test our algorithm on a ubuntu 64bit OS.

The RCC based tracker is adapted from the method developed in <br>
	`"J. S. Stahl and S. Wang, “Edge grouping combining boundary and region information,” IEEE Trans. Image Processing, vol. 16, no. 10, pp. 2590–2606, 2007."` <br>
We download the code from https://cse.sc.edu/~songwang/software.html.

INSTALLATION
-------------------------------------------------
Git clone or Download the zip file and unzip it.<br>
Go to the root path<br>
	`"make"`<br>
then, use command "./SalientBoundaryTracking video_proc gp_flag video_path" to run the code.<br>
For example, read video from .avi file:<br>
	`"./SalientBoundaryTracking 1 1 ./TEST_VIDEO.avi"`<br>
or read video from webcam:<br>
	 `"./SalientBoundaryTracking 0 1"`<br>

video_proc: `"0"`: read file from your webcam, `"1"`: read video from a .avi file <br> 
gp_flag: 	`"0"`: use RRC adapted method to track,  `"1"`: use our BDSP to track <br>
video_path: specify the video sequence to be loaded. If video_proc == 0, it will be neglected.

TRACKING CONTROL
--------------------------------------------------
(1). Click mouse left button to initialize a<br>
　　polygon in `"Display"` window.<br>
　　When read video from a webcam,<br>
　　you have to press key `"i"` before<br>
　　initialization.<br>
(2). Click mouse middle button to finish initialization.<br>
(3). Press key `"space"` to run tracking continuously or<br>
　　press key `"s"` to run tracking frame by frame.<br>
　　Key `"r"` can recover the<br>
　　continuously tracking model).`"Esc"` will quit the tracking.<br>
(4). Tracking output are two or three files:<br>
　　`"TRACKING_RESULTS.txt"`<br>
　　`"RECORDED_RESULTS.avi"`<br>
　　`"RECOREDED_VIDEO.avi"`(if video_proc==0).<br>
