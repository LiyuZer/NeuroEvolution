CC=clang++

evolution: main.cpp Network.hpp Node.hpp organism.hpp Utilities.hpp
	$(CC) --std=c++17 -O3 -o evolution main.cpp


.PHONY: clean
clean:
	rm evolution
