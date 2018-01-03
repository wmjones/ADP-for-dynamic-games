VPATH = src:latex
CXX = g++-7

all : figs/plots.gif figs/%.jpg plot_ANN_Ind_Func.png plot_g.png
	afplay /System/Library/Components/CoreAudio.component/Contents/SharedSupport/SystemSounds/system/payment_success.aif -v .5

figs/plots.gif : figs/%.jpg
	convert -set delay '%[fx:t==0 ? 30 : 40 - t/(n-1)]' -loop 0 figs/%d_plot.jpeg[1-301] plots.gif

figs/%.jpg : model_data.csv plotting.r
	Rscript src/plotting.r

plot_ANN_Ind_Func.png : plotting_ANN_Ind_Func.r
	Rscript src/plotting_ANN_Ind_Func.r

plot_g.png : plotting_g.r
	Rscript src/plotting_g.r

# Cheb_Ind_Func.png

# Chebyshev_Polynomials_of_the_First_Kind

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

pdf : final.tex
	cd latex;\
	ls;\
	pdflatex final;\
	bibtex final;\
	pdflatex final;\
	pdflatex final;\
	cd ../;

pdfclean :
	cd latex;\
	rm *.bbl;\
	rm *.blg;\
	rm *.log;\
	rm *.out;\
	rm *.aux;\
	cd ../;

view : final.pdf
	cd latex;\
	open final.pdf;\
	cd ../;

clean :
	rm program
	rm -rf program.dSYM
	rm figs/*.jpg
	rm model_data.csv
	rm figs/plots.gif
	afplay /System/Library/Components/CoreAudio.component/Contents/SharedSupport/SystemSounds/system/payment_success.aif -v .5
