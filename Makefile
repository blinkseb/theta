#order matters
DIRS = src libconfig liblbfgs plugins root bin test

SVN_TAG := $( a=`svn info | grep URL` )
DUMMY := $( echo svn tag: $(SVN_TAG) )

all:
	@for d in $(DIRS); do $(MAKE) -C $$d || break; done

clean:
	@for d in $(DIRS); do ( [ -d $$d ] && $(MAKE) -C $$d clean ) || true; done
	@$(MAKE) -C test/test-stat clean

#ok, this is not very nice, as it hardcodes the documentation copy target, but as long as I am the only
# developer, it should work well:
doc:
	@doxygen
	@( warn=`wc -l doxygen-warn.txt | cut -f 1 -d" "`; if [ $$warn -gt 0 ]; then echo There have been about $$warn warnings from doxygen, see doxygen-warn.txt; fi )
	@cp doc/tabs.css doc/html
	@SVN_TAG=`svn info | grep URL:`; SVN_TAG=$${SVN_TAG##*/}; echo svn tag: $$SVN_TAG
	@if [ "`hostname`" = "ekplx22" && "$$SVN_TAG" = "stable" ]; then rsync -a --del doc/* /usr/users/ott/public_html/theta; fi
	@if [ "`hostname`" = "ekplx22" && "$$SVN_TAG" = "testing" ]; then rsync -a --del doc/* /usr/users/ott/public_html/theta/testing; fi

.PHONY: clean all doc

