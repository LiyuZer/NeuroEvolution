CC=clang++

evolution: main.cpp Network.hpp Node.hpp organism.hpp Utilities.hpp dna.hpp
	$(CC) --std=c++17 -Wall -O3 -o  evolution main.cpp


.PHONY: clean
clean:
	rm evolution
