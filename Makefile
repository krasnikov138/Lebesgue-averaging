CXX      := g++
CXXFLAGS := --std=c++17
TARGET   := main
OBJ_DIR  := run/bin
APP_DIR  := run
SRC      := $(wildcard src/*.cpp)
OBJECTS  := $(SRC:src/%.cpp=$(OBJ_DIR)/%.o)

all: $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: src/%.cpp
    mkdir -p $(OBJ_DIR)
    $(CXX) $(CXXFLAGS) -c $< -o $@

$(APP_DIR)/$(TARGET): $(OBJECTS)
    $(CXX) $(CXXFLAGS) -o $(APP_DIR)/$(TARGET) $^

debug: CXXFLAGS += -DDEBUG -g
debug: all

release: CXXFLAGS += -O3
release: all

clean:
    rm -rvf $(OBJ_DIR)/*
    rm -vf $(APP_DIR)/$(TARGET)
