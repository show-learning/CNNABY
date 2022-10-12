all:
	g++ L1_to_L4.cpp -O0 -fno-stack-protector -laby -lencrypto_utils -lrelic_s -lotextension -lstdc++fs -lboost_system -std=c++17 -pthread -L/usr/local/lib -lcrypto -lgmp -lgmpxx -lssl -lrt -lSHR -o L12L4
	g++ L5_for.cpp -O0 -fno-stack-protector -laby -lencrypto_utils -lrelic_s -lotextension -lstdc++fs -lboost_system -std=c++17 -pthread -L/usr/local/lib -lcrypto -lgmp -lgmpxx -lssl -lrt -lSHR -o L5_for
	g++ L6_to_L8.cpp -O0 -fno-stack-protector -laby -lencrypto_utils -lrelic_s -lotextension -lstdc++fs -lboost_system -std=c++17 -pthread -L/usr/local/lib -lcrypto -lgmp -lgmpxx -lssl -lrt -lSHR -o L62L8
L_four:
	./L12L4 127.0.0.1 7766 SERVER /home/wei/ABY/bin/circ/ > L4_output.txt
run:
	./L5_for 127.0.0.1 7766 SERVER /home/wei/ABY/bin/circ/ > L5_output_1.txt
	./L5_for 127.0.0.1 7766 SERVER /home/wei/ABY/bin/circ/ > L5_output_2.txt
	./L5_for 127.0.0.1 7766 SERVER /home/wei/ABY/bin/circ/ > L5_output_3.txt
	./L5_for 127.0.0.1 7766 SERVER /home/wei/ABY/bin/circ/ > L5_output_4.txt
	./L5_for 127.0.0.1 7766 SERVER /home/wei/ABY/bin/circ/ > L5_output_5.txt
	./L5_for 127.0.0.1 7766 SERVER /home/wei/ABY/bin/circ/ > L5_output_6.txt
	./L5_for 127.0.0.1 7766 SERVER /home/wei/ABY/bin/circ/ > L5_output_7.txt
	./L5_for 127.0.0.1 7766 SERVER /home/wei/ABY/bin/circ/ > L5_output_8.txt
	./L5_for 127.0.0.1 7766 SERVER /home/wei/ABY/bin/circ/ > L5_output_9.txt
	./L5_for 127.0.0.1 7766 SERVER /home/wei/ABY/bin/circ/ > L5_output_10.txt
	./L5_for 127.0.0.1 7766 SERVER /home/wei/ABY/bin/circ/ > L5_output_11.txt
	./L5_for 127.0.0.1 7766 SERVER /home/wei/ABY/bin/circ/ > L5_output.txt
out:
	./L62L8 127.0.0.1 7766 SERVER /home/wei/ABY/bin/circ/
clean:
	rm -rf L12L4 L5_for L62L8 test.yuv
