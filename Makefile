# Attempt to load a config.make file.
# If none is found, project defaults in config.project.make will be used.
ifneq ($(wildcard config.make),)
	include config.make
endif

# make sure the the OF_ROOT location is defined
ifndef OF_ROOT
	OF_ROOT=../../..
endif

# this project's dir (config.make normally sets it; default if built without config.make)
ifndef PROJECT_PATH
	PROJECT_PATH := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
endif

# fail early with a clear message if openFrameworks isn't where OF_ROOT points
ifeq ($(wildcard $(OF_ROOT)/libs/openFrameworksCompiled/project/makefileCommon/compile.project.mk),)
$(error openFrameworks 0.12.1 not found at OF_ROOT=$(OF_ROOT) — install it, or run: make OF_ROOT=/path/to/of_v0.12.1_osx_release)
endif

# call the project makefile!
include $(OF_ROOT)/libs/openFrameworksCompiled/project/makefileCommon/compile.project.mk

# --- project dependencies (trackplot) -----------------------------------------
# Build the vendored libBigWig static lib and fetch the external ofxDropdown addon
# automatically before the app links. Order-only prerequisites ('| deps') so that,
# once present, they never force a relink. Run them by hand with `make deps`.
LIBBIGWIG   := $(PROJECT_PATH)/libBigWig/libBigWig.a
OFXDROPDOWN := $(OF_ROOT)/addons/ofxDropdown

.PHONY: deps
deps: $(LIBBIGWIG) $(OFXDROPDOWN)

$(LIBBIGWIG):
	$(MAKE) -C $(PROJECT_PATH)/libBigWig libBigWig.a

$(OFXDROPDOWN):
	@test -d "$(OF_ROOT)/addons" || { echo "ERROR: OF_ROOT=$(OF_ROOT) has no addons/ dir — install openFrameworks 0.12.1 (or run: make OF_ROOT=/path/to/of_v0.12.1_osx_release)"; exit 1; }
	git clone --depth 1 https://github.com/roymacdonald/ofxDropdown "$(OFXDROPDOWN)"

Release Debug: | deps
