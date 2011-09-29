include ./Makefile.options

#order matters
DIRS = src libconfig liblbfgs plugins root bin test

ifeq ($(OPTIONS_BUILD_LLVM),yes)
   DIRS += llvm
endif

all:
	@for d in $(DIRS); do ( cd $$d; $(MAKE) ) || break; done

clean:
	@for d in $(DIRS); do ( cd $$d && $(MAKE) clean ) || true; done
	@cd test/test-stat; $(MAKE) clean

#ok, this is not very nice, as it hardcodes the documentation copy target, but as long as I am the only
# developer, it should work well:
doc:
	@doxygen
	@( warn=`wc -l doxygen-warn.txt | cut -f 1 -d" "`; if [ $$warn -gt 0 ]; then echo There have been about $$warn warnings from doxygen, see doxygen-warn.txt; fi )
	@cp doc/tabs.css doc/html
	@if [ "`hostname`" = "ekplx22" ]; then \
	    SVN_TAG=`svn info | grep URL:`; SVN_TAG=$${SVN_TAG##*/}; echo svn tag: $$SVN_TAG; \
	   if [ "$$SVN_TAG" = "stable" ]; then rsync -a --del doc/* /usr/users/ott/public_html/theta; fi;\
	   if [ "$$SVN_TAG" = "testing" ]; then rsync -a --del doc/* /usr/users/ott/public_html/theta/testing; fi;\
	fi
   

.PHONY: clean all doc

