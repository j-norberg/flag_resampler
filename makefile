
# same for all
CXX := g++

# Incl. directories
FLAGS := -Wall -Wshadow -O3 -march=corei7-avx
TARGET := flag_resampler

$(TARGET) : main.o pffft_double.o streamer.o writer.o file_reader.o 
	$(CXX) -o $(TARGET) main.o pffft-double.o streamer.o writer.o file_reader.o

main.o: src/main.cpp src/file_reader.h src/writer.h src/streamer.h
	$(CXX) $(FLAGS) -c src/main.cpp

streamer.o: src/streamer.cpp src/streamer.h src/sample_producer.h
	$(CXX) $(FLAGS) -c src/streamer.cpp

writer.o: src/writer.cpp src/writer.h src/sample_producer.h
	$(CXX) $(FLAGS) -c src/writer.cpp

file_reader.o: src/file_reader.cpp src/file_reader.h src/sample_producer.h
	$(CXX) $(FLAGS) -c src/file_reader.cpp

pffft_double.o: src/dep/pffft-double.c src/dep/pffft-double.h
	$(CXX) $(FLAGS) -c src/dep/pffft-double.c

clean:
	-@rm -vf $(TARGET)
	-@rm -rvf *.o
