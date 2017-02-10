EXECS= salida.x
MPICC?=mpicc

placas.pdf : data.txt grafica.py
	python grafica.py

data.txt : salida.x
	./salida.x > data.txt

salida.x : placas.c
	gcc placas.c -o salida.x

clean:
	rm -f salida.x data.txt
