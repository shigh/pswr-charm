
CXX = g++
TLDLIBS += -lboost_unit_test_framework

test: testsuite
	./testsuite --log_level=test_suite

testsuite: $(TOBJS) $(OBJS)
	$(CXX) -o $@ $^ $(LDLIBS) $(TLDLIBS)

$(OBJS): %.o: ../%.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@

$(TOBJS): %.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@

# For emacs flymake mode
check-syntax:
	$(CXX) $(CXXFLAGS) -c -o nul -Wall $(CHK_SOURCES)


