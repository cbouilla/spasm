build: build/Makefile
	$(MAKE) --no-print-directory -C build

build/Makefile:
	mkdir build ; cd build ; cmake .. ; cd ..

mrproper:
	rm -rf build

clean:
	$(MAKE) --no-print-directory -C build clean

check: build/Makefile
	$(MAKE) --no-print-directory -C build check

longtests: build/Makefile
	$(MAKE) --no-print-directory -C build longtests

.PHONY: build mrproper clean check longtests
