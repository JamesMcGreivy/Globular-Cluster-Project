CC=g++


#GENERAL COMMAND
plots : 
	@make test1plots
	@make test2plots
	@make test3plots
	@make test4plots

#All make commands required for test case 1

#Compiles the integrators using the test1.h file
./test1/boris.out : ./integrators/boris.cpp main.cpp ./test1/test1.h
	@$(CC) -w -o ./test1/boris.out main.cpp -include ./integrators/boris.cpp -include ./test1/test1.h

./test1/collocation.out : ./integrators/collocation.cpp main.cpp ./test1/test1.h
	@$(CC) -w -o ./test1/collocation.out main.cpp -include ./integrators/collocation.cpp -include ./test1/test1.h

#Runs the simulations
./test1/data/boris100.csv : ./test1/boris.out
	@./test1/boris.out 0.01 ./test1/data/boris100.csv

./test1/data/boris1000.csv : ./test1/boris.out
	@./test1/boris.out 0.001 ./test1/data/boris1000.csv

./test1/data/collocation2.csv : ./test1/collocation.out
	@./test1/collocation.out 0.5 ./test1/data/collocation2.csv

#Makes the plots
./test1/energy_errors.eps : ./test1/graph_energy.py ./test1/data/boris100.csv ./test1/data/boris1000.csv ./test1/data/collocation2.csv
	@python3 ./test1/graph_energy.py 

./test1/phase_errors.eps : ./test1/graph_phase_error.py ./test1/data/boris100.csv ./test1/data/boris1000.csv ./test1/data/collocation2.csv
	@python3 ./test1/graph_phase_error.py 

#General commands
test1data :
	@make ./test1/data/boris100.csv
	@make ./test1/data/boris1000.csv
	@make ./test1/data/collocation2.csv
test1plots : 
	@make ./test1/energy_errors.eps
	@make ./test1/phase_errors.eps



#All make commands required for test case 2

#Compiles the integrators using the test2.h file
./test2/boris.out : ./integrators/boris.cpp main.cpp ./test2/test2.h
	@$(CC) -w -o ./test2/boris.out main.cpp -include ./integrators/boris.cpp -include ./test2/test2.h

./test2/collocation.out : ./integrators/collocation.cpp main.cpp ./test2/test2.h
	@$(CC) -w -o ./test2/collocation.out main.cpp -include ./integrators/collocation.cpp -include ./test2/test2.h

#Runs the simulations
./test2/data/boris100.csv : ./test2/boris.out
	@./test2/boris.out 0.01 ./test2/data/boris100.csv

./test2/data/boris1000.csv : ./test2/boris.out
	@./test2/boris.out 0.001 ./test2/data/boris1000.csv

./test2/data/collocation5.csv : ./test2/collocation.out
	@./test2/collocation.out 0.2 ./test2/data/collocation5.csv

./test2/data/collocation10.csv : ./test2/collocation.out
	@./test2/collocation.out 0.1 ./test2/data/collocation10.csv

#Makes the plots
./test2/energy_errors.eps : ./test2/graph_energy.py ./test2/data/boris100.csv ./test2/data/boris1000.csv ./test2/data/collocation5.csv ./test2/data/collocation10.csv
	@python3 ./test2/graph_energy.py

./test2/x_trajectory.eps : ./test2/graph_x.py ./test2/data/boris100.csv ./test2/data/boris1000.csv ./test2/data/collocation5.csv ./test2/data/collocation10.csv
	@python3 ./test2/graph_x.py

#General Commands
test2data :
	@make ./test2/data/boris100.csv
	@make ./test2/data/boris1000.csv
	@make ./test2/data/collocation5.csv
	@make ./test2/data/collocation10.csv

test2plots :
	@make ./test2/energy_errors.eps
	@make ./test2/x_trajectory.eps



#All make commands required for test case 3

#Compiles the integrators using the test3.h file
./test3/boris.out : ./integrators/boris.cpp main.cpp ./test3/test3.h
	@$(CC) -w -o ./test3/boris.out main.cpp -include ./integrators/boris.cpp -include ./test3/test3.h

./test3/collocation.out : ./integrators/collocation.cpp main.cpp ./test3/test3.h
	@$(CC) -w -o ./test3/collocation.out main.cpp -include ./integrators/collocation.cpp -include ./test3/test3.h

#Runs the simulations
./test3/data/boris10.csv : ./test3/boris.out
	@./test3/boris.out 0.1 ./test3/data/boris10.csv

./test3/data/boris100.csv : ./test3/boris.out
	@./test3/boris.out 0.01 ./test3/data/boris100.csv

./test3/data/boris1000.csv : ./test3/boris.out
	@./test3/boris.out 0.001 ./test3/data/boris1000.csv

./test3/data/collocation2.csv : ./test3/collocation.out
	@./test3/collocation.out 0.5 ./test3/data/collocation2.csv

./test3/data/collocation5.csv : ./test3/collocation.out
	@./test3/collocation.out 0.2 ./test3/data/collocation5.csv

./test3/data/collocation10.csv : ./test3/collocation.out
	@./test3/collocation.out 0.1 ./test3/data/collocation10.csv

#Makes the plots
./test3/energy_errors.eps : ./test3/graph_energy.py ./test3/data/boris10.csv ./test3/data/boris100.csv ./test3/data/boris1000.csv ./test3/data/collocation2.csv ./test3/data/collocation5.csv ./test3/data/collocation10.csv
	@python3 ./test3/graph_energy.py

./test3/mag_moment.eps : ./test3/graph_magnetic_moment.py ./test3/data/boris10.csv ./test3/data/boris100.csv ./test3/data/boris1000.csv ./test3/data/collocation2.csv ./test3/data/collocation5.csv ./test3/data/collocation10.csv
	@python3 ./test3/graph_magnetic_moment.py

./test3/trajectory.eps : ./test3/graph_trajectory.py ./test3/data/boris10.csv ./test3/data/boris100.csv ./test3/data/boris1000.csv ./test3/data/collocation2.csv ./test3/data/collocation5.csv ./test3/data/collocation10.csv
	@python3 ./test3/graph_trajectory.py

./test3/trajectory_error.eps : ./test3/trajectory_error.py ./test3/data/boris10.csv ./test3/data/boris100.csv ./test3/data/boris1000.csv ./test3/data/collocation2.csv ./test3/data/collocation5.csv ./test3/data/collocation10.csv
	@python3 ./test3/trajectory_error.py

#General Commands
test3data :
	@make ./test3/data/boris10.csv
	@make ./test3/data/boris100.csv
	@make ./test3/data/boris1000.csv
	@make ./test3/data/collocation2.csv
	@make ./test3/data/collocation5.csv
	@make ./test3/data/collocation10.csv

test3plots :
	@make ./test3/energy_errors.eps
	@make ./test3/mag_moment.eps
	@make ./test3/trajectory.eps
	@make ./test3/trajectory_error.eps


#All make commands required for test case 4

#Compiles the integrators using the test4.h file and different values for epsilon
./test4/eps8.out : ./integrators/collocation.cpp main.cpp ./test4/test4.h ./test4/eps/eps8.h
	@$(CC) -w -o ./test4/eps8.out main.cpp -include ./test4/eps/eps8.h -include ./integrators/collocation.cpp -include ./test4/test4.h 

./test4/eps7.out : ./integrators/collocation.cpp main.cpp ./test4/test4.h ./test4/eps/eps7.h
	@$(CC) -w -o ./test4/eps7.out main.cpp -include ./test4/eps/eps7.h -include ./integrators/collocation.cpp -include ./test4/test4.h 

./test4/eps6.out : ./integrators/collocation.cpp main.cpp ./test4/test4.h ./test4/eps/eps6.h
	@$(CC) -w -o ./test4/eps6.out main.cpp -include ./test4/eps/eps6.h -include ./integrators/collocation.cpp -include ./test4/test4.h 

./test4/eps5.out : ./integrators/collocation.cpp main.cpp ./test4/test4.h ./test4/eps/eps5.h
	@$(CC) -w -o ./test4/eps5.out main.cpp -include ./test4/eps/eps5.h -include ./integrators/collocation.cpp -include ./test4/test4.h 

./test4/eps4.out : ./integrators/collocation.cpp main.cpp ./test4/test4.h ./test4/eps/eps4.h
	@$(CC) -w -o ./test4/eps4.out main.cpp -include ./test4/eps/eps4.h -include ./integrators/collocation.cpp -include ./test4/test4.h 

./test4/eps3.out : ./integrators/collocation.cpp main.cpp ./test4/test4.h ./test4/eps/eps3.h
	@$(CC) -w -o ./test4/eps3.out main.cpp -include ./test4/eps/eps3.h -include ./integrators/collocation.cpp -include ./test4/test4.h 

./test4/eps2.out : ./integrators/collocation.cpp main.cpp ./test4/test4.h ./test4/eps/eps2.h
	@$(CC) -w -o ./test4/eps2.out main.cpp -include ./test4/eps/eps2.h -include ./integrators/collocation.cpp -include ./test4/test4.h 

./test4/eps1.out : ./integrators/collocation.cpp main.cpp ./test4/test4.h ./test4/eps/eps1.h
	@$(CC) -w -o ./test4/eps1.out main.cpp -include ./test4/eps/eps1.h -include ./integrators/collocation.cpp -include ./test4/test4.h 

./test4/eps0.out : ./integrators/collocation.cpp main.cpp ./test4/test4.h ./test4/eps/eps0.h
	@$(CC) -w -o ./test4/eps0.out main.cpp -include ./test4/eps/eps0.h -include ./integrators/collocation.cpp -include ./test4/test4.h 


#Runs the simulations
#T = 1
./test4/data/1eps8.csv : ./test4/eps8.out
	@./test4/eps8.out 1 ./test4/data/1eps8.csv

./test4/data/1eps7.csv : ./test4/eps7.out
	@./test4/eps7.out 1 ./test4/data/1eps7.csv

./test4/data/1eps6.csv : ./test4/eps6.out
	@./test4/eps6.out 1 ./test4/data/1eps6.csv

./test4/data/1eps5.csv : ./test4/eps5.out
	@./test4/eps5.out 1 ./test4/data/1eps5.csv

./test4/data/1eps4.csv : ./test4/eps4.out
	@./test4/eps4.out 1 ./test4/data/1eps4.csv

./test4/data/1eps3.csv : ./test4/eps3.out
	@./test4/eps3.out 1 ./test4/data/1eps3.csv

./test4/data/1eps2.csv : ./test4/eps2.out
	@./test4/eps2.out 1 ./test4/data/1eps2.csv

./test4/data/1eps1.csv : ./test4/eps1.out
	@./test4/eps1.out 1 ./test4/data/1eps1.csv

#T = 1/2
./test4/data/2eps8.csv : ./test4/eps8.out
	@./test4/eps8.out 0.5 ./test4/data/2eps8.csv

./test4/data/2eps7.csv : ./test4/eps7.out
	@./test4/eps7.out 0.5 ./test4/data/2eps7.csv

./test4/data/2eps6.csv : ./test4/eps6.out
	@./test4/eps6.out 0.5 ./test4/data/2eps6.csv

./test4/data/2eps5.csv : ./test4/eps5.out
	@./test4/eps5.out 0.5 ./test4/data/2eps5.csv

./test4/data/2eps4.csv : ./test4/eps4.out
	@./test4/eps4.out 0.5 ./test4/data/2eps4.csv

./test4/data/2eps3.csv : ./test4/eps3.out
	@./test4/eps3.out 0.5 ./test4/data/2eps3.csv

./test4/data/2eps2.csv : ./test4/eps2.out
	@./test4/eps2.out 0.5 ./test4/data/2eps2.csv

./test4/data/2eps1.csv : ./test4/eps1.out
	@./test4/eps1.out 0.5 ./test4/data/2eps1.csv

#T = 1/5
./test4/data/5eps8.csv : ./test4/eps8.out
	@./test4/eps8.out 0.2 ./test4/data/5eps8.csv

./test4/data/5eps7.csv : ./test4/eps7.out
	@./test4/eps7.out 0.2 ./test4/data/5eps7.csv

./test4/data/5eps6.csv : ./test4/eps6.out
	@./test4/eps6.out 0.2 ./test4/data/5eps6.csv

./test4/data/5eps5.csv : ./test4/eps5.out
	@./test4/eps5.out 0.2 ./test4/data/5eps5.csv

./test4/data/5eps4.csv : ./test4/eps4.out
	@./test4/eps4.out 0.2 ./test4/data/5eps4.csv

./test4/data/5eps3.csv : ./test4/eps3.out
	@./test4/eps3.out 0.2 ./test4/data/5eps3.csv

./test4/data/5eps2.csv : ./test4/eps2.out
	@./test4/eps2.out 0.2 ./test4/data/5eps2.csv

./test4/data/5eps1.csv : ./test4/eps1.out
	@./test4/eps1.out 0.2 ./test4/data/5eps1.csv


#Makes the plots
./test4/e_vs_eps.eps : ./test4/plot_eps_vs_e.py ./test4/data/5eps8.csv ./test4/data/5eps7.csv ./test4/data/5eps6.csv ./test4/data/5eps5.csv ./test4/data/5eps4.csv ./test4/data/5eps3.csv ./test4/data/5eps2.csv ./test4/data/5eps1.csv ./test4/data/1eps8.csv ./test4/data/1eps7.csv ./test4/data/1eps6.csv ./test4/data/1eps5.csv ./test4/data/1eps4.csv ./test4/data/1eps3.csv ./test4/data/1eps2.csv ./test4/data/1eps1.csv ./test4/data/2eps8.csv ./test4/data/2eps7.csv ./test4/data/2eps6.csv ./test4/data/2eps5.csv ./test4/data/2eps4.csv ./test4/data/2eps3.csv ./test4/data/2eps2.csv ./test4/data/2eps1.csv
	@python3 ./test4/plot_eps_vs_e.py


#General commands
test4plots : 
	@make ./test4/e_vs_eps.eps