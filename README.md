# VFA-for-dynamic-games
To get this project working on your own computer running OSX first open terminal and run
```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```
then run
```
binaries=(
    gcc
    cmake
    openblas --with-openmp
    nlopt
    dlib --cc=gcc-7 --with-openblas
)
brew install ${binaries[@]}
```
Then simply cd to the directory where the makefile is on your computer and run 
```make all```
