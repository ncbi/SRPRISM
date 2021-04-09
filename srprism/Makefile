# UNCOMMENT THE FOLLOWING LINE TO ENABLE DIRECT SRA ACCESS
# export USE_SRA = 1

# UNCOMMENT THE FOLLOWING LINES TO POINT TO EXISTING NGS INSTALLATION
# AND FILL IN THE CORRECT PATHS
#
# export NGS_HOME = 
# export VDB_HOME =

ifdef USE_SRA

ifndef NGS_HOME

NGS_DIR := $(CURDIR)/NGS
NGSDEP = $(NGS_DIR)/ngs.done

VDB_PATH := $(NGS_DIR)/vdb_out
NGS_PATH := $(NGS_DIR)/ngs_out
BUILD_PATH := $(NGS_DIR)/build

VDB_INCL := -I $(VDB_PATH)/include
VDB_LIB := -L $(VDB_PATH)/lib64
#VDB_LIB := -L $(VDB_PATH)/lib64 -Xlinker -rpath=$(VDB_PATH)/lib64
NGS_INCL := -I $(NGS_PATH)/include
NGS_LIB := -L $(NGS_PATH)/lib64

else

VDB_INCL := -I $(VDB_HOME)/include
VDB_LIB := -L $(VDB_HOME)/lib
NGS_INCL := -I $(NGS_HOME)/include
NGS_LIB := -L $(NGS_HOME)/lib

endif

endif

export INC = $(NGS_INCL) $(VDB_INCL)
export LIB = $(NGS_LIB) $(VDB_LIB)

SUBDIRS = lib app

all: $(NGSDEP)
	for i in $(SUBDIRS) ; do $(MAKE) -C $$i ; done

ifdef USE_SRA

ifndef NGS_HOME

$(NGSDEP):
	rm -fr $(NGS_DIR) ;
	mkdir -p $(NGS_DIR)/ngs ;
	mkdir $(BUILD_PATH) ;
	mkdir $(NGS_PATH) ;
	mkdir $(VDB_PATH) ;
	cd $(NGS_DIR)/ngs; git init; git remote add -f origin https://github.com/ncbi/ngs.git; git config core.sparseCheckout true; echo "ngs-sdk" >> .git/info/sparse-checkout; git pull origin master ;
	cd $(NGS_DIR)/ngs/ngs-sdk; ./configure --prefix=$(NGS_PATH) --build-prefix=$(BUILD_PATH); make; make install ;
	cd $(NGS_DIR); git clone https://github.com/ncbi/ncbi-vdb.git ;
	cd $(NGS_DIR)/ncbi-vdb; ./configure --prefix=$(VDB_PATH) --build-prefix=$(BUILD_PATH); make; make install ;
	touch $@

endif

endif

clean:
	for i in $(SUBDIRS) ; do cd $$i && $(MAKE) clean ; cd .. ; done

