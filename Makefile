# definitions
objRRC = SalientBoundaryTracking.o construct_graph.o construct_graph_dt.o  mat_main.o match.o CmdLine.o BDSP_optimization.o CommonFunctions.o shapeTracker.o
srcRRC = SalientBoundaryTracking.cpp construct_graph.cpp construct_graph_dt.cpp mat_main.cpp match.c CmdLine.cpp BDSP_optimization.cpp CommonFunctions.cpp shapeTracker.cpp

#linker to use
lnk = g++
#compiler to use
cc = gcc
#uncomment for debugging
dbg = -g -Wall

# MAKE it happen

all: SalientBoundaryTracking

SalientBoundaryTracking: $(objRRC)
	$(lnk) $(dbg) -o SalientBoundaryTracking $(objRRC) EDLinesLib.a `pkg-config --libs opencv`

$(objRRC): $(srcRRC)
	$(cc) $(dbg) `pkg-config --cflags opencv` -c $(srcRRC)

clean:
	@rm -f $(objRRC) SalientBoundaryTracking


#all:
#	gcc `pkg-config --cflags opencv` -o VideoEDLines $(srcRRC) EDLinesLib.a libopencv_imgproc.so.2.4.5 libopencv_core.so.2.4.5 `pkg-config --libs opencv`

