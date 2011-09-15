include ../Makefile.options

libsources := $(wildcard *.cpp)
libobjects := $(patsubst %.cpp,.bin/%.o,$(libsources))
deps := $(patsubst %.cpp,.bin/%.d,$(libsources))

ifeq ($(OPTIONS_USE_CRLIBM),yes)
   libobjects += .bin/log.o
   deps += .bin/log.d
endif

libobjects += .bin/log2_dot.o

include ../Makefile.rules

all: ../lib/libtheta.so



../lib/libtheta.so: $(libobjects)
	@$(LINK_SO) -o $@ $+

.bin/%.o: ../crlibm/%.c
	@$(COMPILE) -x c -o $@ $<


# switch between implementations of log2_dot
ifeq ($(OPTIONS_GENERIC_ARCH),yes)
.bin/log2_dot.o: log2_dot.cxx
	@$(COMPILE) -o $@ $<
else
.bin/log2_dot.o: log2_dot.s
	@$(COMPILE) -o $@ $<
endif

clean-subdir:
	@rm -f ../lib/libtheta.so
