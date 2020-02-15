.SECONDARY:

MAKEFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
MAKEFILE_DIR := $(dir $(MAKEFILE_PATH))

# Pseudotargets that control compilation
IS_DEBUG = $(findstring debug,$(MAKECMDGOALS))
IS_UNOPTIMIZED = $(findstring unoptimized,$(MAKECMDGOALS))

# Install dir
PREFIX = /usr/local
INSTALL_BIN = $(PREFIX)/bin

# Compiler: try clang++, fall back to g++
CPP = clang++
ifeq (, $(shell which $(CPP)))
CPP = g++
endif

# Boost
BOOST_PROGRAM_OPTIONS = program_options
BOOST_OBJ_FILES =
# Try to figure out where Boost is
# NB pkg-config support for Boost is lacking; see https://svn.boost.org/trac/boost/ticket/1094
BOOST_PREFIX = /usr
ifeq (,$(wildcard $(BOOST_PREFIX)/include/boost/program_options.hpp))
BOOST_PREFIX = /usr/local
ifeq (,$(wildcard $(BOOST_PREFIX)/include/boost/program_options.hpp))
BOOST_PREFIX =
endif
endif
BOOST_FLAGS =
BOOST_LIBS =
ifneq (,$(BOOST_PREFIX))
BOOST_FLAGS := -I$(BOOST_PREFIX)/include
BOOST_LIBS := -L$(BOOST_PREFIX)/lib -lboost_$(BOOST_PROGRAM_OPTIONS)
endif

# Compiler & linker flags
ALL_FLAGS = $(BOOST_FLAGS)
ALL_LIBS = $(BOOST_LIBS)

ifneq (,$(IS_DEBUG))
CPP_FLAGS = -std=c++11 -g -DUSE_VECTOR_GUARDS -DDEBUG
else
ifneq (,$(IS_UNOPTIMIZED))
CPP_FLAGS = -std=c++11 -g
else
CPP_FLAGS = -std=c++11 -g -O3
endif
endif
CPP_FLAGS += $(ALL_FLAGS) -Isrc
LD_FLAGS = -lstdc++ -lm $(ALL_LIBS)

# Files
CPP_FILES = $(wildcard src/*.cpp)
OBJ_FILES = $(subst src/,obj/,$(subst .cpp,.o,$(CPP_FILES)))

# pwd
PWD = $(shell pwd)

# /bin/sh
SH = /bin/sh

# Targets
CARNAVAL = carnaval

all: $(CARNAVAL)

install: $(CARNAVAL)
	cp bin/$(CARNAVAL) $(INSTALL_BIN)/$(CARNAVAL)

# Main build rules
bin/%: $(OBJ_FILES) obj/%.o
	@test -e $(dir $@) || mkdir -p $(dir $@)
	$(CPP) $(LD_FLAGS) -o $@ obj/$*.o $(OBJ_FILES)

obj/%.o: src/%.cpp
	@test -e $(dir $@) || mkdir -p $(dir $@)
	$(CPP) $(CPP_FLAGS) -c -o $@ $<

obj/%.o: target/%.cpp
	@test -e $(dir $@) || mkdir -p $(dir $@)
	$(CPP) $(CPP_FLAGS) -c -o $@ $<

$(CARNAVAL): bin/$(CARNAVAL)

clean:
	rm -rf bin/$(CARNAVAL) obj/*

# Fake pseudotargets
debug unoptimized:
