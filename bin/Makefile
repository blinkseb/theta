include ../Makefile.options

deps := .bin/merge.d .bin/theta.d

include ../Makefile.rules

all: merge theta


merge: .bin/merge.o
	@$(LINK_EXE) -o $@ $<

theta: .bin/theta.o
	@$(LINK_EXE) -o $@ $<

clean-subdir:
	@rm -f merge theta

