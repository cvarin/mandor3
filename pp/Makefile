COMPILER=icc
CFLAGS=-std=c99 -c -O3 -I$$HOME/bin/fftw3/include
LFLAGS=-lfftw3 -lm -L$$HOME/bin/fftw3/lib
OUTPUT=pp.out

all: compile link

compile:
	$(COMPILER) $(CFLAGS) source/cfgreader.c -o o/cfgreader.o
	$(COMPILER) $(CFLAGS) source/data_loader.c -o o/data_loader.o
	$(COMPILER) $(CFLAGS) source/main.c -o o/main.o
	$(COMPILER) $(CFLAGS) source/matrix.c -o o/matrix.o
	$(COMPILER) $(CFLAGS) source/plot.c -o o/plot.o
	$(COMPILER) $(CFLAGS) source/processor.c -o o/processor.o
	$(COMPILER) $(CFLAGS) source/report.c -o o/report.o
	$(COMPILER) $(CFLAGS) source/utils.c -o o/utils.o
	$(COMPILER) $(CFLAGS) source/io_utils.c -o o/io_utils.o

link:
	$(COMPILER) $(LFLAGS) o/*.o $$HOME/bin/fftw3/lib/libfftw3.a -o $(OUTPUT)


clean:
	rm -f o/*.*
	rm -f dat/*.*
	rm -f pics/*.*
	rm -f plot.plt
	rm -f *.stackdump
	rm -f $(OUTPUT)
