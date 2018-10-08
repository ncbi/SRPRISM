HEADERS =   cfasta_stream.hpp \
            cfastq_stream.hpp \
            fasta_stream.hpp \
            fastq_stream.hpp \
            paired_stream.hpp \
            seqdef.hpp \
            seqinput.hpp \
            seqinput_factory.hpp \
            seqinput_multistream.hpp \
            seqinput_sam.hpp \
            seqinput_sra.hpp \
            stream_base.hpp \
            stream_factory.hpp \

SOURCES =   cfasta_stream.cpp \
            cfastq_stream.cpp \
            fasta_stream.cpp \
            fastq_stream.cpp \
            seqdef.cpp \
            seqinput_factory.cpp \
            seqinput_multistream.cpp \
            seqinput_sam.cpp \
            seqinput_sra.cpp \
            stream_factory.cpp \

OBJECTS = $(SOURCES:.cpp=.o)

INCLUDES = -I.. -I../..
CFLAGS = -std=c++11 -O3 -Wno-deprecated-declarations

all: libseq.a

libseq.a : $(OBJECTS) $(HEADERS)
	ar cr $@ $(OBJECTS) ;
	ranlib $@

%.o: %.cpp $(HEADERS)
	g++ $(INCLUDES) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.a *.o


