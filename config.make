# openFrameworks 0.12.1 install. Not bundled in the repo — install it separately.
# Override per machine:  make OF_ROOT=/path/to/of_v0.12.1_osx_release
OF_ROOT ?= /Applications/of_v0.12.1_osx_release

# This project's dir, resolved from config.make's own location (works from any clone path).
PROJECT_PATH := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))

# --- libBigWig (built locally with -DNOCURL: local .bw files only) ---
PROJECT_CFLAGS  = -I$(PROJECT_PATH)/libBigWig
PROJECT_LDFLAGS = -L$(PROJECT_PATH)/libBigWig -lBigWig -lz
