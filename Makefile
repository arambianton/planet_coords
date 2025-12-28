CXX      := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -pedantic
DBGFLAGS := -std=c++17 -g   -Wall -Wextra -pedantic
TARGET   := solar
SRC      := solar.cpp
OBJ      := $(SRC:.cpp=.o)

.PHONY: all debug run clean

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

debug: CXXFLAGS = $(DBGFLAGS)
debug: clean $(TARGET)

run: $(TARGET)
	./$(TARGET) $(ARGS)

clean:
	rm -f $(OBJ) $(TARGET)
