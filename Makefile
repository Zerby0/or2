tsp: *.c
	gcc *.c -lm -o tsp

clean:
	rm -f tsp
