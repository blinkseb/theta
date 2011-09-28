libsources := $(wildcard *.cpp)
libobjects := $(patsubst %.cpp,.bin/%.o,$(libsources))
deps := $(patsubst %.cpp,.bin/%.d,$(libsources))

include ../Makefile.options
include ../Makefile.rules

EXTRA_COMPILE += $(shell llvm-config --cppflags)
EXTRA_LINK_SO += $(shell llvm-config --ldflags) -lLLVM-2.7

all: ../lib/llvm-plugins.so

../lib/llvm-plugins.so: $(libobjects)
	@$(LINK_SO) -o $@ $+

clean-subdir:
	@rm -f ../lib/llvm-plugins.so
