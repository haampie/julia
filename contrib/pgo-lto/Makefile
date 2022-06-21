.PHONY: top clean

STAGE0_BUILD:=$(CURDIR)/stage0.build
STAGE1_BUILD:=$(CURDIR)/stage1.build
STAGE2_BUILD:=$(CURDIR)/stage2.build

STAGE0_TOOLS:=$(STAGE0_BUILD)/usr/tools/

PROFILE_DIR:=$(CURDIR)/profiles
PROFILE_FILE:=$(PROFILE_DIR)/merged.prof
JULIA_ROOT:=$(CURDIR)/../..
CXXFILT:=c++filt

LLVM_PROFDATA:=$(STAGE0_TOOLS)llvm-profdata
export CC:=$(STAGE0_TOOLS)clang
export CXX:=$(STAGE0_TOOLS)clang++
export LD:=$(STAGE0_TOOLS)ld.lld
export AR:=$(STAGE0_TOOLS)llvm-ar
export RANLIB:=$(STAGE0_TOOLS)llvm-ranlib

# When building a single libLLVM.so we need to increase -vp-counters-per-site
# significantly
COUNTERS_PER_SITE:=6

AFTER_STAGE1_MESSAGE:=Consider removing $(PROFILE_DIR) to start with a clean slate. $\
    Then run Julia to collect realistic profile data, for example: \`julia -O3 -e $\
    'using Pkg; Pkg.add("Unitful"); Pkg.test("Unitful")'\`. This should produce about $\
	15MB of data in $(PROFILE_DIR). Note that running extensive scripts may result in $\
	counter overflows, which can be detected by running \`make $(PROFILE_FILE)\`. $\
	Afterwards run `make stage2`.

$(STAGE0_BUILD) $(STAGE1_BUILD) $(STAGE2_BUILD):
	$(MAKE) -C $(JULIA_ROOT) O=$@ configure

stage0: export CFLAGS=
stage0: export CXXFLAGS=
stage0: export LDFLAGS=
stage0: export USE_BINARYBUILDER_LLVM=1
stage0: | $(STAGE0_BUILD)
	$(MAKE) -C $(STAGE0_BUILD)/deps install-clang install-llvm install-llvm-tools
	# Turn [cd]tors into init/fini_array sections in libclang_rt, since lld
	# doesn't do that, and otherwise the profile constructor is not executed
	find $< -name 'libclang_rt.profile-*.a' -exec objcopy --rename-section .ctors=.init_array --rename-section .dtors=.fini_array {} +
	touch $@

$(STAGE1_BUILD): stage0
stage1: export CFLAGS=-fprofile-generate=$(PROFILE_DIR) -Xclang -mllvm -Xclang -vp-counters-per-site=$(COUNTERS_PER_SITE)
stage1: export CXXFLAGS=-fprofile-generate=$(PROFILE_DIR) -Xclang -mllvm -Xclang -vp-counters-per-site=$(COUNTERS_PER_SITE)
stage1: export LDFLAGS=-fuse-ld=lld -flto=thin -fprofile-generate=$(PROFILE_DIR)
stage1: export USE_BINARY_BUILDER=0
stage1: | $(STAGE1_BUILD)
	$(MAKE) -C $(STAGE1_BUILD) && echo $(AFTER_STAGE1_MESSAGE) && touch $@

stage2: export CFLAGS=-fprofile-use=$(PROFILE_FILE)
stage2: export CXXFLAGS=-fprofile-use=$(PROFILE_FILE)
stage2: export LDFLAGS=-fuse-ld=lld -flto=thin -fprofile-use=$(PROFILE_FILE) -Wl,--icf=safe
stage2: export USE_BINARY_BUILDER=0
stage2: $(PROFILE_FILE) | $(STAGE2_BUILD)
	$(MAKE) -C $(STAGE2_BUILD) && touch $@

$(PROFILE_FILE): stage1
	$(LLVM_PROFDATA) merge -output=$@ $(wildcard $(PROFILE_DIR)/*.profraw)

# show top 50 functions
top: $(PROFILE_FILE)
	$(LLVM_PROFDATA) show --topn=50 $< | $(CXXFILT)

clean:
	rm -f stage0 stage1 stage2 $(PROFILE_FILE)
