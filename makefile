CC=g++

make:main.cpp
	@$(CC) -w main.cpp $(epsilon)
	@./a.out
clean:
	@rm a.out
