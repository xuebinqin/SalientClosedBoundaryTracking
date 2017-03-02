Slient Closed Boundary Tracking
============================
This code is implemeted for the papper `"Slaient Closed Boundary Tracking via Line Segments Perceptual Grouping, Xuebin Qin, Shida He, Camilo Perez Quintero, Abhineet Singh, masood dehghan and Martin Jagersand."` which is submited to IROS 2017.

This implementation are tested on `Opencv 2.4.9(3.1.0)`, `Boost 1.63.0` and `ubuntu 14.04 64 bit`.

The line detector used here is EDLines(EDLines.a) which is proposed in `"C. Akinlar and C. Topal, EDLines: A real-time line segment detector with a false detection control, Pattern Recognition Letters, vol. 32, no. 13, pp. 1633-1642, 2011."` and can be downloaded from http://ceng.anadolu.edu.tr/CV/downloads/downloads.aspx. We include the lib in the root file. It is worth to note that we are using the 64 bit ubuntu version of EDLines. We sugggest to test our algorithm on a ubuntu 64bit OS.

The RCC based tracker is adapted from the method developed in `"J. S. Stahl and S. Wang, “Edge grouping combining boundary and region information,” IEEE Trans. Image Processing, vol. 16, no. 10, pp. 2590–2606, 2007."` We download the code from https://cse.sc.edu/~songwang/software.html.

INSTALLATION
-------------------------------------------------
Git clone or Download the zip file and unzip it.<br>
Go to the root path and "make" it,<br>
then, use command "./SalientBoundaryTracking video_proc gp_flag video_path" to run the code.<br>
`"./SalientBoundaryTracking 1 1 ./TEST_VIDEO.avi"`<br>

video_proc: <br>
            "0": read file from your webcam
            "1": read video from a .avi file
gp_flag: 	<br>
             "0": use RRC adapted method to track
             "1": use our BDSP to track
video_path: specify the video sequence to be loaded.<br>
            If video_proc == 0, it will be neglected.<br>

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
     "TRACKING_RESULTS.txt"<br>
     "RECORDED_RESULTS.avi"<br>
     "RECOREDED_VIDEO.avi"(if video_proc==0).<br>

# SlientClosedBoundaryTracking
