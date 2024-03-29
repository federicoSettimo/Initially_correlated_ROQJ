roqj.o: roqj.h roqj.cpp
	g++ roqj.cpp -c -o roqj.o -std=c++20 -O3 -ffast-math -fno-math-errno

qubit_ph_cov: roqj.o Examples/qubit_ph_cov.cpp
	g++ Examples/qubit_ph_cov.cpp roqj.o -o Examples/qubit_ph_cov.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/qubit_ph_cov.x

dephasing_nM: Examples/dephasing_nM.cpp
	g++ Examples/dephasing_nM.cpp -o Examples/dephasing_nM.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/dephasing_nM.x
	python3 Examples/plot_deph_nM.py

dephasing_nM_partial_info: Examples/dephasing_nM_partial_info.cpp
	g++ Examples/dephasing_nM_partial_info.cpp -o Examples/dephasing_nM_partial_info.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/dephasing_nM_partial_info.x
	python3 Examples/plot_deph_nM.py

error: Examples/dephasing_nM_error.cpp Examples/dephasing_nM_error_partial_info.cpp
	g++ Examples/dephasing_nM_error.cpp -o Examples/dephasing_nM_error.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/dephasing_nM_error.x
	g++ Examples/dephasing_nM_error_partial_info.cpp -o Examples/dephasing_nM_error_partial_info.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/dephasing_nM_error_partial_info.x
	python3 Examples/plot_deph_nM_error.py

dephasing_2qubits: Examples/dephasing_2qubits.cpp
	g++ Examples/dephasing_2qubits.cpp -o Examples/dephasing_2qubits.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/dephasing_2qubits.x
	python3 Examples/plot_deph_2qubits.py

Jaynes-Cummings_env_n: Examples/Jaynes-Cummings_env_n.cpp
	g++ Examples/Jaynes-Cummings_env_n.cpp -o Examples/Jaynes-Cummings_env_n.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/Jaynes-Cummings_env_n.x
	python3 Examples/plot_JC.py