HEADERS = options_parser.hpp options_parser_priv.hpp
SOURCES = options_parser.cpp srprism.cpp

OBJECTS = $(SOURCES:.cpp=.o)
LIBRARIES = ../lib/srprism/libsrprism.a \
            ../lib/seq/libseq.a \
            ../lib/common/libcommon.a

INCLUDES = -I. -I.. -I../lib
#CFLAGS = -std=c++11 -g -O0 -Wno-deprecated-declarations
CFLAGS = -std=c++11 -DNDEBUG -O3 -Wno-deprecated-declarations

ifdef USE_SRA

#LIBS = -pthread -lncbi-ngs-c++ -lngs-c++ -lncbi-vdb-static -ldl 
#LIBS = -pthread -lncbi-ngs-c++ -lngs-c++ -lncbi-vdb -ldl 
LIBS = -pthread -lncbi-ngs-c++ -lngs-c++ -lncbi-ngs-static -lncbi-vdb-static -ldl 

endif

LIBS += -lpthread -lz -lbz2

all: srprism

srprism: $(OBJECTS) $(HEADERS) $(LIBRARIES)
	g++ $(LIB) -o $@ $(OBJECTS) $(LIBRARIES) $(LIBS)

%.o: %.cpp $(HEADERS)
	g++ $(INCLUDES) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.a *.o srprism
