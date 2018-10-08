SUBDIRS = common seq srprism

all:
	for i in $(SUBDIRS) ; do $(MAKE) -C $$i ; done

clean:
	for i in $(SUBDIRS) ; do cd $$i && $(MAKE) clean ; cd .. ; done

