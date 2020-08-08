CC=g++

make:
	@$(CC) -w -o run.out  main.cpp
run:
	@./run.out $(config)
clean:
	@rm run.out
