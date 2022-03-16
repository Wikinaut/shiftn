# shiftn
ShiftN - Python

This repo aims towards a Python command-line version of ShiftN https://shiftn.de originally written by Marcus Hebel. The full source code is available on his website.


#### Step 1: Line detection (burns.cpp, burns.py)

How to compile burns.cpp (see [INSTALL.txt](https://github.com/Wikinaut/shiftn/edit/main/INSTALL.txt))

```
sudo apt install g++
sudo apt install cmake
sudo apt install libopencv-dev
 
cmake .
make
./burns

Das erzeugt dann aus Ullsteinhaus.jps → Ullsteinhaus_Linien.png und Ullsteinhaus_Linien.txt, und man sieht zumindest, dass der Code funktioniert.
```
