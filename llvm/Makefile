libsources := $(wildcard *.cpp)
libobjects := $(patsubst %.cpp,.bin/%.o,$(libsources))
deps := $(patsubst %.cpp,.bin/%.d,$(libsources))

include ../Makefile.options
include ../Makefile.rules

EXTRA_COMPILE += $(shell llvm-config --cppflags)

all: ../lib/llvm-plugins.so

../lib/llvm-plugins.so: $(libobjects)
	@$(LINK_SO) $(shell llvm-config --libfiles nativecodegen engine scalaropts ipo)

clean-subdir:
	@rm -f ../lib/llvm-plugins.so
