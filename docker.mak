deregnet: lemon-headers grbfrc 
	docker build -t sebwink/deregnet .

grbfrc: gurobi-cpp 
	cd grbfrc && docker build -t sebwink/grbfrc .

gurobi-cpp:
	cd docker/gurobi-cpp && make

lemon-headers:
	cd docker/lemon-headers && make
