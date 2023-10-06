TARGET_EXEC1 = compareFingerPrints
TARGET_EXEC2 = compareToLibrary
TARGET_EXEC3 = getFingerPrint

SRCS1 = src/point.cpp src/writhe.cpp src/mainFileCompareFingerprints.cpp
SRCS2 = src/point.cpp src/writhe.cpp src/mainFileCompareFingerPrintAgainstLibrary.cpp
SRCS3 = src/point.cpp src/writhe.cpp src/mainFileFingerprint.cpp

OBJS1 = $(SRCS1:.cpp=.o)
OBJS2 = $(SRCS2:.cpp=.o)
OBJS3 = $(SRCS3:.cpp=.o)

CXX = g++
CXXFLAGS = -O3 -std=gnu++17
RM = rm -f


all: $(TARGET_EXEC1) $(TARGET_EXEC2) $(TARGET_EXEC3)

$(TARGET_EXEC1): $(OBJS1)
	$(CXX) $(CXXFLAGS) -o $(TARGET_EXEC1) $(OBJS1) 

$(TARGET_EXEC2): $(OBJS2)
	$(CXX) $(CXXFLAGS) -o $(TARGET_EXEC2) $(OBJS2) 

$(TARGET_EXEC3): $(OBJS3)
	$(CXX) $(CXXFLAGS) -o $(TARGET_EXEC3) $(OBJS3) 

source/cpp/.cpp.o:
	$(CXX) $(CXXFLAGS) $< -o $@ 

clean:
	$(RM) $(OBJS2) $(OBJS1) $(OBJS3) $(TARGET_EXEC1) $(TARGET_EXEC2) $(TARGET_EXEC3)

