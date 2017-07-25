all : plots.gif
	afplay /System/Library/Components/CoreAudio.component/Contents/SharedSupport/SystemSounds/system/payment_success.aif -v .5

plots.gif : %.jpg
	convert -set delay '%[fx:t==0 ? 30 : 40 - t/(n-1)]' -loop 0 plot%d.jpg[1-100] plots.gif

%.jpg : model_data.csv		# not sure why this runs if all images already created
	Rscript plotting.r

model_data.csv : program
	./program

program : main.cpp parameters.cpp model.cpp fitting.cpp
	# g++-7 -g -Wall -fopenmp -lstdc++ -lm -std=c++11 -lnlopt -I/usr/local/include/ main.cpp model.cpp parameters.cpp fitting.cpp -o program
	g++-7 -O3 -fopenmp -lstdc++ -lm -std=c++11 -lnlopt -ldlib -I/usr/local/include/ main.cpp model.cpp parameters.cpp fitting.cpp -o program

# build :
# 	g++-7 main.cpp model.cpp parameters.cpp fitting.cpp -o program -std=c++11 -lmlpack -larmadillo -lboost_serialization -lboost_program_options -lnlopt -fopenmp -I/usr/local/include/
# 	afplay /System/Library/Components/CoreAudio.component/Contents/SharedSupport/SystemSounds/system/payment_success.aif -v .5
# -march=native

run : program
	./program
	afplay /System/Library/Components/CoreAudio.component/Contents/SharedSupport/SystemSounds/system/payment_success.aif -v .5

clean :
	rm program
	rm *.jpg
	rm model_data.csv
	rm plots.gif
	afplay /System/Library/Components/CoreAudio.component/Contents/SharedSupport/SystemSounds/system/payment_success.aif -v .5
