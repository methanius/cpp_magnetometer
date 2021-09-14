SHELL = bash
CXXFLAGS = -Wall -O3 -std=c++17
LDLIBS = -larmadillo -lm 
CC = g++-7

OBJECTS = main.o ehrenfest_chain.o make_operators.o choi_jamiolkowski.o homodyne_emission.o homodyne_PQS.o
SOURCES = main.cpp ehrenfest_chain.cpp make_operators.cpp choi_jamiolkowski.cpp homodyne_emission.cpp homodyne_PQS.cpp


.PHONY:default clean profile

default:out.txt
	cat $<

out.txt:main variables.txt example_images beta_images Makefile
#	./$< variables.txt > $@
#	python3 cpp_OU_figures.py 

beta_images:main beta_images_variables.txt
	index=0; 
	for b in {0..4}; \
		do \
		for bb in {1..9..2}; \
			do \
			sed -i "16s/.*/$$b.$$bb/g" beta_images_variables.txt; \
				folder_name="beta_images/beta_$$b.$$bb" \
				index=$$((index+1)); \
				./$< beta_images_variables.txt $$folder_name `sed "$${index}q;d" seeds_beta_figures.txt`; \
			done \
		done
	#python3 function_of_beta_curve.py

example_images:main example_images_variables.txt
	#The integer at the end here is the seed
	./$< example_images_variables.txt $@ 201505292
	python3 cpp_OU_figures.py $@

main: $(OBJECTS)


clean:
	$(RM) out.txt main $(OBJECTS)

profile:
	for i in $(SOURCES); do \
		$(CC) $(CPPFLAGS) -pg -c $$i; \
	done
	
	$(CC) $(OBJECTS) $(LDLIBS) -pg -o main
	./main
	gprof main > profile.stats


