#order matters
DIRS = src libconfig liblbfgs plugins root bin test

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
	@if [ "`hostname`" = "ekplx22" ]; then rsync -a --del doc/* /usr/users/ott/public_html/theta; fi

.PHONY: clean all doc

