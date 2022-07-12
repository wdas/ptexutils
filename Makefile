uname_S ?= $(shell uname)
uname_R ?= $(shell uname -r | sed -e s,-.*,,)
uname_M ?= $(shell uname -m)
FLAVOR ?= optimize

prefix_base ?= $(uname_S)-$(uname_R)-$(uname_M)-$(FLAVOR)
prefix ?= $(CURDIR)/$(prefix_base)
CMAKE_BUILD_DIR ?= build/$(prefix_base)

# Add more variables here to force a rebuild when they change
TRACK_FLAGS = $(prefix):$(FLAVOR)

all: install

# Record build flags so that we can force cmake to reconfigure
$(CMAKE_BUILD_DIR)/FLAGS: FORCE
	mkdir -p $(CMAKE_BUILD_DIR)
	@FLAGS='$(TRACK_FLAGS)'; \
	if test x"$$FLAGS" != x"`cat $(CMAKE_BUILD_DIR)/FLAGS 2>/dev/null`" ; then \
		echo >&2 "-- New build flags"; \
		echo "$$FLAGS" > $(CMAKE_BUILD_DIR)/FLAGS; \
	fi

clean:
	rm -rf -- $(CMAKE_BUILD_DIR) Linux-*

$(CMAKE_BUILD_DIR)/Makefile: CMakeLists.txt src/CMakeLists.txt $(CMAKE_BUILD_DIR)/FLAGS
	cd $(CMAKE_BUILD_DIR) && cmake -DCMAKE_INSTALL_PREFIX=$(prefix) ../..

cmake: $(CMAKE_BUILD_DIR)/Makefile

install: cmake
	$(MAKE) --no-print-directory -C $(CMAKE_BUILD_DIR) install $(MFLAGS)

.PHONY: all install cmake
.PHONY: FORCE
