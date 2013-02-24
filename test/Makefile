include ../Makefile.options

ifeq ($(OPTIONS_BUILD_TESTS),yes)

# build test*.cpp as binary test:
testsources := $(wildcard test*.cpp utils.cpp)
testobjects := $(patsubst %.cpp,.bin/%.o,$(testsources))

# build test*.cxx as .so; will be loaded as plugin by tests:
libtestsources := $(wildcard test*.cxx)
libtestobjects := $(patsubst %.cxx,.bin/%.o,$(libtestsources))

deps := $(patsubst %.o,%.d,$(libtestobjects) $(testobjects))

include ../Makefile.rules

all: ../bin/test ../lib/liblibtest.so

../bin/test: $(testobjects)
	@$(LINK_EXE)

../lib/liblibtest.so: $(libtestobjects)
	@$(LINK_SO)

else

all:
	@echo -e \\nNOTE: not building tests as this was switched off in Makefile.options\\n

deps:=
include ../Makefile.rules

endif

clean-subdir:
	@rm -f ../bin/test ../lib/liblibtest.so

