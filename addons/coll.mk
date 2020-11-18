CC	= g++
FF	= g77
LD	= g++


BINARIES = colcount
HSOURCES = 
HEADERS  = $(HSOURCES:.cxx=.h)
SOURCES  = $(HSOURCES) colcount.cxx
OBJECTS  = $(SOURCES:.cxx=.o) 

CPPOPT = -Wno-deprecated `root-config --cflags`
COPTN = -c -O3 ${CPPOPT} 
COPT = -c -O3 ${CPPOPT}  
LLIB = `root-config --libs` -L$(ROOTSYS)/lib -lm  -lgcc
LOPT = 

all: $(BINARIES)

colcount: $(OBJECTS)
	$(LD) $(LOPT) -o $@ $(OBJECTS) $(LLIB)

colcount.o: colcount.cxx
	$(CC) $(COPT) -o $@ $^

%.o: %.cxx
	$(CC) -o $@ $< $(COPT)

