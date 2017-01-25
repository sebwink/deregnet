.PHONY: all drgnt avgdrgnt clean destroy

all : drgnt avgdrgnt build

build :
	mkdir -f build

drgnt :
	make -f drgnt.mak

avgdrgnt : grbfrc
	#make -f avgdrgnt.mak

grbfrc :
	cd grbfrc
	make 
	cd ..

clean :
	rm -f build/*.o

destroy :
	rm -f build/*.o
	rm -f bin/drgnt bin/avgdrgnt
