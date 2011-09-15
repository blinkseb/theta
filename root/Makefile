include ../Makefile.options

ifeq ($(shell which root-config || echo -n A),A)
all:
	@echo -e \\nNOTE: no root-config found. Not building root-dependent stuff.\\n

include ../Makefile.rules

else

sources := $(wildcard *.cpp)
objects := $(patsubst %.cpp,.bin/%.o,$(sources))
deps := $(patsubst %.cpp,.bin/%.d,$(sources)) .bin/histos.d .bin/db2root.d

include ../Makefile.rules

# just use -Wall here, do not warn about old-stryle casts, etc. Otherwise, we will be swamped by warnings
# from root.
WFLAGSCPP := -Wall
EXTRA_COMPILE := $(shell root-config --cflags)
EXTRA_LINK_SO := $(shell root-config --libs) -lMinuit2 $(shell root-config --ldflags)
EXTRA_LINK_EXE := $(EXTRA_LINK_SO)

all: ../lib/root.so ../bin/histos

ifeq ($(OPTIONS_BUILD_TESTS),yes)
all: create_testhistos create_testtree
endif

../bin/histos: .bin/histos.o
	@$(LINK_EXE) -o $@ $<

create_testhistos: .bin/create_testhistos.o
	@$(LINK_EXE) -o $@ $<

create_testhistos2: .bin/create_testhistos2.o
	@$(LINK_EXE) -o $@ $<

create_histos: .bin/create_histos.o
	@$(LINK_EXE) -o $@ $<

create_testtree: .bin/create_testtree.o
	@$(LINK_EXE) -o $@ $<

../lib/root.so: $(objects)
	@$(LINK_SO) -o $@ $+

endif

#define clean target in any case:
clean-subdir:
	@rm -f histos create_testhistos ../lib/root.so
