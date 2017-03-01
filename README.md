This code is implemeted for the papper "Slaient Closed Boundary Tracking via Line Segments Perceptual Grouping, Xuebin Qin, Shida He, Camilo Perez Quintero, Abhineet Singh and Martin Jagersand." which is submited to IROS 2017.

This implementation are tested on Opencv 2.4.9(3.1.0), Boost 1.63.0 and ubuntu 14.04 64 bit.

The line detector used here is EDLines(EDLines.a) which is proposed in "C. Akinlar and C. Topal, EDLines: A real-time line segment detector with a false detection control, Pattern Recognition Letters, vol. 32, no. 13, pp. 1633-1642, 2011." and can be downloaded from http://ceng.anadolu.edu.tr/CV/downloads/downloads.aspx. We include the lib in the root file. It is worth to note that we are using the 64 bit ubuntu version of EDLines. We sugggest to test our algorithm on a ubuntu 64bit OS.

The RCC based tracker is adapted from the method developed in "J. S. Stahl and S. Wang, “Edge grouping combining boundary and region information,” IEEE Trans. Image Processing, vol. 16, no. 10, pp. 2590–2606, 2007." We download the code from https://cse.sc.edu/~songwang/software.html.


Git clone or Download the zip file and unzip it.
Go to the root path and "make" it,
then, use command "./SalientBoundaryTracking <video_proc> <gp_flag> <video_path>" to run the code, 
such as "./SalientBoundaryTracking 1 1 ./TEST_VIDEO.avi".
-------------------------------------------------
video_proc: 
            "0": read file from your webcam
            "1": read video from a .avi file
gp_flag: 
             "0": use RRC adapted method to track
             "1": use our BDSP to track
video_path: specify the video sequence to be loaded.
            If video_proc == 0, it will be neglected.
-------------------------------------------------
(1). Click mouse left button to initialize a
     polygon in "Display" window.
     When read video from a webcam,
     you have to press key "i" before
     initialization.
(2). Click mouse middle button to finish initialization.
(3). Press key "space" to run tracking continuously or
     press key "s" to run tracking frame by frame.
     Key "r" can recover the
     continuously tracking model).
(4). Tracking output are two or three files: 
     "TRACKING_RESULTS.txt"
     "RECORDED_RESULTS.avi"
     "RECOREDED_VIDEO.avi"(if video_proc==0).
-------------------------------------------------

# SlientClosedBoundaryTracking
