.PHONY: start generate clean clean_object_file execute

all: start generate

object_file = Q2.out
code = Q2.cpp

start: 
	@echo "We are running the make file.... \n"

generate: ${code}
	@echo "Compiling the code.... \n"
	g++ -std=c++11 ${code} -o ${object_file}

execute: ${object_file}
	@echo "Executing the code: \n"
	./${object_file}

clean:
	@echo "Cleaning up the generated text files.... \n"
	rm *.csv
	
clean_object_file:
	@echo "Removing the object file.... \n"
	rm ${object_file}
