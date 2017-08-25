VPATH = src:latex
CXX = g++-7

all : figs/plots.gif
	afplay /System/Library/Components/CoreAudio.component/Contents/SharedSupport/SystemSounds/system/payment_success.aif -v .5

figs/plots.gif : figs/%.jpg
	convert -set delay '%[fx:t==0 ? 30 : 40 - t/(n-1)]' -loop 0 figs/%d_plot.jpeg[1-301] plots.gif

figs/%.jpg : model_data.csv		# not sure why this runs if all images already created
	Rscript src/plotting.r

model_data.csv : program
	./program

program : main.cpp parameters.cpp model.cpp fitting.cpp
	$(CXX) -march=native -O3 -w -fopenmp -lstdc++ -lm -std=c++11 -lnlopt -ldlib -I/usr/local/include/ $^ -o $@

d : main.cpp parameters.cpp model.cpp fitting.cpp
	$(CXX) -g -Wall -fopenmp -lstdc++ -lm -std=c++11 -lnlopt -ldlib -I/usr/local/include/ $^ -o $@
	./$@
	rm $@
	rm -rf $@.dSYM
	afplay /System/Library/Components/CoreAudio.component/Contents/SharedSupport/SystemSounds/system/payment_success.aif -v .5

run : program
	./program
	afplay /System/Library/Components/CoreAudio.component/Contents/SharedSupport/SystemSounds/system/payment_success.aif -v .5

plot : model_data.csv
	Rscript src/plotting.r

gif :
	convert -set delay '%[fx:t==0 ? 30 : 40 - t/(n-1)]' -loop 0 figs/%d_plot.jpeg[1-301] plots.gif

pdf : draft1.tex
	cd latex;\
	ls;\
	pdflatex draft1;\
	bibtex draft1;\
	pdflatex draft1;\
	pdflatex draft1;\
	cd ../;

view : draft1.pdf
	cd latex;\
	open draft1.pdf;\
	cd ../;

clean :
	rm program
	rm -rf program.dSYM
	rm figs/*.jpg
	rm model_data.csv
	rm figs/plots.gif
	afplay /System/Library/Components/CoreAudio.component/Contents/SharedSupport/SystemSounds/system/payment_success.aif -v .5
