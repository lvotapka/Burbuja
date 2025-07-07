#!/bin/bash

for bubble_name in */; do
	cd "$bubble_name"
	bubble_name="${bubble_name%/}"
	python ../../find_bubbles.py "$bubble_name".pdb > cpu_time
	python ../../find_bubbles_cuda.py "$bubble_name".pdb > cuda_time
	#mv bubbles.pdb bubbles_0.7_cut.pdb
	cd ..
done
