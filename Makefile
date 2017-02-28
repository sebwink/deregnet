.PHONY: clean destroy

all :
	./build.sh

clean :
	rm -f build/*.o

destroy :
	rm -f build/*.o
	rm -f bin/drgnt bin/avgdrgnt
	rm -f grbfrc/lib/*
