ifneq ($(MAKECMDGOALS),configure)
	include config.mk
endif

configure:src/configure.ac
	@echo "  AUTOCONF   $@"
	@autoconf --warnings=error -o $@ $< && rm -r autom4te.cache 2>/dev/null



ReadTomographyRAW: src/ReadTomographyRAW.cpp
	$(CPP) -O3 $(OPT) $(LD_OPT) $^ -o $@
	 