.PHONY: top clean

STAGE1_BUILD:=$(CURDIR)/stage1.build
STAGE2_BUILD:=$(CURDIR)/stage2.build

PROFILE_DIR:=$(CURDIR)/profiles
PROFILE_FILE:=$(PROFILE_DIR)/merged.prof
JULIA_ROOT:=$(CURDIR)/../..
CXXFILT:=c++filt

# Note: can't use BB's LLVM right now; -fprofile-generate compiles, but doesn't
# generate profiles at runtime
LLVM_PROFDATA:=llvm-profdata
export CC:=clang
export CXX:=clang++
export LD:=ld.lld
export AR:=llvm-ar
export RANLIB:=llvm-ranlib
export USECLANG:=1
export USE_BINARYBUILDER_LLVM:=0
# export VERBOSE:=1

# When building a single libLLVM.so we need to increase -vp-counters-per-site
# significantly
COUNTERS_PER_SITE:=6

AFTER_STAGE1_MESSAGE:=Consider removing $(PROFILE_DIR) to start with a clean slate. $\
    Then run Julia to collect realistic profile data, for example: \`julia -O3 -e $\
    'using Pkg; Pkg.add("Unitful"); Pkg.test("Unitful")'\`. This should produce about $\
	15MB of data in $(PROFILE_DIR). Note that running extensive scripts may result in $\
	counter overflows, which can be detected by running \`make $(PROFILE_FILE)\`. $\
	Afterwards run `make stage2`.

$(STAGE1_BUILD) $(STAGE2_BUILD):
	$(MAKE) -C $(JULIA_ROOT) O=$@ configure

stage1: export CFLAGS=-fprofile-generate=$(PROFILE_DIR) -Xclang -mllvm -Xclang -vp-counters-per-site=$(COUNTERS_PER_SITE)
stage1: export CXXFLAGS=-fprofile-generate=$(PROFILE_DIR) -Xclang -mllvm -Xclang -vp-counters-per-site=$(COUNTERS_PER_SITE)
stage1: export LDFLAGS=-fuse-ld=lld -flto=thin -fprofile-generate=$(PROFILE_DIR)
stage1: | $(STAGE1_BUILD)
	$(MAKE) -C $(STAGE1_BUILD) && echo Consider removing $(PROFILE_DIR) to start with a clean slate && touch $@

stage2: export CFLAGS=-fprofile-use=$(PROFILE_FILE)
stage2: export CXXFLAGS=-fprofile-use=$(PROFILE_FILE)
stage2: export LDFLAGS=-fuse-ld=lld -flto=thin -fprofile-use=$(PROFILE_FILE)
stage2: $(PROFILE_FILE) | $(STAGE2_BUILD)
	$(MAKE) -C $(STAGE2_BUILD) && touch $@

$(PROFILE_FILE): stage1
	$(LLVM_PROFDATA) merge -output=$@ $(wildcard $(PROFILE_DIR)/*.profraw)

# show top 50 functions
top: $(PROFILE_FILE)
	$(LLVM_PROFDATA) show --topn=50 $< | $(CXXFILT)

clean:
	rm -f stage1 stage2 $(PROFILE_FILE)
