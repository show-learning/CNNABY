all:
	g++ L4.cpp -O0 -fno-stack-protector -laby -lencrypto_utils -lrelic_s -lotextension -lstdc++fs -lboost_system -std=c++17 -pthread -L/usr/local/lib -lcrypto -lgmp -lgmpxx -lssl -lrt -lSHR -o L4
	g++ simfilter.cpp -O0 -fno-stack-protector -laby -lencrypto_utils -lrelic_s -lotextension -lstdc++fs -lboost_system -std=c++17 -pthread -L/usr/local/lib -lcrypto -lgmp -lgmpxx -lssl -lrt -lSHR -o imfltr
L_four:
	./L4 127.0.0.1 7766 SERVER /home/wei/ABY/bin/circ/ > L4_output.txt
run:
	time ./imfltr 127.0.0.1 7766 SERVER /home/wei/ABY/bin/circ/ > L5_output.txt

clean:
	rm -rf L4 imfltr test.yuv
