libsources := $(wildcard *.cpp)
libobjects := $(patsubst %.cpp,.bin/%.o,$(libsources))
deps := $(patsubst %.cpp,.bin/%.d,$(libsources))

include ../Makefile.options
include ../Makefile.rules

all: ../lib/core-plugins.so

../lib/core-plugins.so: $(libobjects)
	@$(LINK_SO) -o $@ $+

clean-subdir:
	@rm -f ../lib/core-plugins.so
