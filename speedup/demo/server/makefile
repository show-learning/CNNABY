all:
	g++ S_share.cpp -fno-stack-protector -laby -lencrypto_utils -lrelic_s -lotextension -lstdc++fs -lboost_system -std=c++17 -pthread -L/usr/local/lib -lcrypto -lgmp -lgmpxx -lssl -lrt -lSHR -o S_share
	g++ S_value.cpp -fno-stack-protector -laby -lencrypto_utils -lrelic_s -lotextension -lstdc++fs -lboost_system -std=c++17 -pthread -L/usr/local/lib -lcrypto -lgmp -lgmpxx -lssl -lrt -lSHR -o S_value
	g++ quick.cpp -o quick
	g++ R_share.cpp -fno-stack-protector -laby -lencrypto_utils -lrelic_s -lotextension -lstdc++fs -lboost_system -std=c++17 -pthread -L/usr/local/lib -lcrypto -lgmp -lgmpxx -lssl -lrt -lSHR -o R_share
	g++ R_value.cpp -fno-stack-protector -laby -lencrypto_utils -lrelic_s -lotextension -lstdc++fs -lboost_system -std=c++17 -pthread -L/usr/local/lib -lcrypto -lgmp -lgmpxx -lssl -lrt -lSHR -o R_value
run:
	./S_share 127.0.0.1 7766 SERVER /home/shower/ABY/bin/circ/ > L1_output.txt
	./S_value 127.0.0.1 7766 SERVER /home/shower/ABY/bin/circ/ 
	./quick input.txt test.txt
	./R_share 127.0.0.1 7766 SERVER /home/shower/ABY/bin/circ/ > Output.txt
	./R_value 127.0.0.1 7766 SERVER /home/shower/ABY/bin/circ/ 

clean:
	rm -rf input.txt L1_output.txt test.yuv Output.txt output.txt input.txt test.txt
	rm -rf S_share S_value quick R_share R_value 
