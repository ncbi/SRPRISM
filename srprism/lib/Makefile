SUBDIRS = common seq srprism

all:
	for i in $(SUBDIRS) ; do cd $$i && $(MAKE) ; cd .. ; done

clean:
	for i in $(SUBDIRS) ; do cd $$i && $(MAKE) clean ; cd .. ; done

