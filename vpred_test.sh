
make sim-safe
cp sim-safe benchmarks/
cd benchmarks/
./sim-safe cc1.alpha -O 1stmt.i
#./sim-safe compress95.alpha words < compress95.in > OUT
#./sim-safe anagram.alpha words < anagram.in > OUT
cd ..