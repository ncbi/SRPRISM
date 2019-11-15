HEADERS =   binfile.hpp \
            bits.hpp \
            bzipfile.hpp \
            exception.hpp \
            file.hpp \
            memqsort.hpp \
            text_formatter.hpp \
            textfile.hpp \
            tmpstore.hpp \
            trace.hpp \
            util.hpp \
            zipfile.hpp \

SOURCES =   binfile.cpp \
            bzipfile.cpp \
            file.cpp \
            text_formatter.cpp \
            textfile.cpp \
            tmpstore.cpp \
            trace.cpp \
            util.cpp \
            zipfile.cpp \

OBJECTS = $(SOURCES:.cpp=.o)

INCLUDES = -I.. -I../..
#CFLAGS = -std=c++11 -g -O0 -Wno-deprecated-declarations
CFLAGS = -std=c++11 -O3 -Wno-deprecated-declarations

all: libcommon.a

libcommon.a : $(OBJECTS) $(HEADERS)
	ar cr $@ $(OBJECTS) ;
	ranlib $@

%.o: %.cpp $(HEADERS)
	g++ $(INCLUDES) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.a *.o

