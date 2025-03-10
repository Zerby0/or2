CFLAGS=-std=c11 -D_POSIX_C_SOURCE=200809L -lm

tsp: *.c
	gcc $(CFLAGS) $^ -o $@

tsp-debug: *.c
	gcc -g $(CFLAGS) $^ -o $@

tsp-release: *.c
	gcc -O3 $(CFLAGS) $^ -o $@

clean:
	rm -f tsp tsp-debug tsp-release
