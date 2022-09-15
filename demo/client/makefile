all:
	g++ L4.cpp -O0 -fno-stack-protector -laby -lencrypto_utils -lrelic_s -lotextension -lstdc++fs -lboost_system -std=c++17 -pthread -L/usr/local/lib -lcrypto -lgmp -lgmpxx -lssl -lrt -lSHR -o L4
	g++ simfilter.cpp -O0 -fno-stack-protector -laby -lencrypto_utils -lrelic_s -lotextension -lstdc++fs -lboost_system -std=c++17 -pthread -L/usr/local/lib -lcrypto -lgmp -lgmpxx -lssl -lrt -lSHR -o imfltr
	g++ L6_to_L8.cpp -O0 -fno-stack-protector -laby -lencrypto_utils -lrelic_s -lotextension -lstdc++fs -lboost_system -std=c++17 -pthread -L/usr/local/lib -lcrypto -lgmp -lgmpxx -lssl -lrt -lSHR -o L62L8
L_four:
	./L4 127.0.0.1 7766 CLIENT /home/wei/ABY/bin/circ/ > L4_output.txt
run:
	./imfltr 127.0.0.1 7766 CLIENT /home/wei/ABY/bin/circ/ > L5_output_1.txt
	./imfltr 127.0.0.1 7766 CLIENT /home/wei/ABY/bin/circ/ > L5_output_2.txt
	./imfltr 127.0.0.1 7766 CLIENT /home/wei/ABY/bin/circ/ > L5_output_3.txt
	./imfltr 127.0.0.1 7766 CLIENT /home/wei/ABY/bin/circ/ > L5_output_4.txt
	./imfltr 127.0.0.1 7766 CLIENT /home/wei/ABY/bin/circ/ > L5_output_5.txt
	./imfltr 127.0.0.1 7766 CLIENT /home/wei/ABY/bin/circ/ > L5_output_6.txt
	./imfltr 127.0.0.1 7766 CLIENT /home/wei/ABY/bin/circ/ > L5_output_7.txt
	./imfltr 127.0.0.1 7766 CLIENT /home/wei/ABY/bin/circ/ > L5_output_8.txt
	./imfltr 127.0.0.1 7766 CLIENT /home/wei/ABY/bin/circ/ > L5_output_9.txt
	./imfltr 127.0.0.1 7766 CLIENT /home/wei/ABY/bin/circ/ > L5_output_10.txt
	./imfltr 127.0.0.1 7766 CLIENT /home/wei/ABY/bin/circ/ > L5_output_11.txt
	./imfltr 127.0.0.1 7766 CLIENT /home/wei/ABY/bin/circ/ > L5_output_12.txt
out:
	./L62L8 127.0.0.1 7766 SERVER /home/wei/ABY/bin/circ/
clean:
	rm -rf L4 imfltr test.yuv