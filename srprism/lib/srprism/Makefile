HEADERS =   align.hpp \
            batch.hpp \
            batch_priv.hpp \
            bnf.hpp \
            idx_reader.hpp \
            idxmap_reader.hpp \
            index_base.hpp \
            index_iterator.hpp \
            inplace_align.hpp \
            memmgr.hpp \
            mkidx.hpp \
            mkidx_pass.hpp \
            nmer_iterator.hpp \
            out_base.hpp \
            out_sam.hpp \
            out_tabular.hpp \
            query_acct.hpp \
            query_data.hpp \
            query_store.hpp \
            query_store_priv.hpp \
            result.hpp \
            rmap.hpp \
            scoring.hpp \
            scratch.hpp \
            search.hpp \
            search_mode.hpp \
            search_pass.hpp \
            search_pass_priv.hpp \
            seqiter.hpp \
            seqstore.hpp \
            seqstore_base.hpp \
            seqstore_factory.hpp \
            seqstore_factory_priv.hpp \
            sidmap.hpp \
            srprismdef.hpp \
            stat.hpp \
            tmpres_mgr.hpp \

SOURCES =   align.cpp \
            batch.cpp \
            idx_reader.cpp \
            idxmap_reader.cpp \
            index_iterator.cpp \
            memmgr.cpp \
            mkidx.cpp \
            mkidx_pass.cpp \
            nmer_iterator.cpp \
            out_sam.cpp \
            out_tabular.cpp \
            query_data.cpp \
            query_store.cpp \
            rmap.cpp \
            search.cpp \
            search_pass.cpp \
            seqstore.cpp \
            seqstore_factory.cpp \
            sidmap.cpp \
            tmpres_mgr.cpp \

OBJECTS = $(SOURCES:.cpp=.o)

INCLUDES = -I.. -I../..
#CFLAGS = -std=c++11 -g -O0 -Wno-deprecated-declarations
CFLAGS = -std=c++11 -DNDEBUG -O3 -Wno-deprecated-declarations

all: libsrprism.a

libsrprism.a : $(OBJECTS) $(HEADERS)
	ar cr $@ $(OBJECTS) ;
	ranlib $@

%.o: %.cpp $(HEADERS)
	g++ $(INCLUDES) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.a *.o


