roqj.o: roqj.h roqj.cpp
	g++ roqj.cpp -c -o roqj.o -std=c++20 -O3 -ffast-math -fno-math-errno

roqj_pop.o: roqj_pop.cpp roqj_pop.h
	g++ roqj_pop.cpp -c -o roqj_pop.o -std=c++20 -O3 -ffast-math -fno-math-errno

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

dephasing_0_discord: Examples/dephasing_0_discord.cpp
	g++ Examples/dephasing_0_discord.cpp -o Examples/dephasing_0_discord.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/dephasing_0_discord.x

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

JC_max_ent: Examples/JC_max_ent.cpp
	g++ Examples/JC_max_ent.cpp -o Examples/JC_max_ent.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/JC_max_ent.x
	python3 Examples/plot_JC_max_ent.py

Ising: Examples/Ising.cpp
	g++ Examples/Ising.cpp -o Examples/Ising.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/Ising.x
	python3 Examples/Ising.py

ph_cov_0_discord: roqj.o roqj_pop.o Examples/ph_cov_0_discord.cpp
	g++ Examples/ph_cov_0_discord.cpp roqj.o roqj_pop.o -o Examples/ph_cov_0_discord.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/ph_cov_0_discord.x
	#python3 Examples/plot_ph_cov.py "Jaynes-Cummings, $ \Phi^0_t(Q_x) $" "$ \operatorname{tr}[X \sigma_x] $"

JC_APO: roqj.o roqj_pop.o Examples/JC_APO.cpp
	g++ Examples/JC_APO.cpp roqj.o roqj_pop.o -o Examples/JC_APO.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/JC_APO.x
	#python3 Examples/plot_JC_APO.py "Jaynes-Cummings with APO, $ \Phi^x_t(Q_x) $, $ \vert \Psi_0\rangle = (\vert0,4\rangle + \vert1,0\rangle)/\sqrt{2} $" "$ \operatorname{tr}[X \sigma_x] $"

JC_APO_continuum: roqj.o roqj_pop.o Examples/JC_APO_continuum.cpp
	g++ Examples/JC_APO_continuum.cpp roqj.o roqj_pop.o -o Examples/JC_APO_continuum.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/JC_APO_continuum.x
	#python3 Examples/plot_JC_APO.py "JC with APO, continuum limit, $ \Phi^x_t(Q_x) $, maximally entangled initial state" "$ \operatorname{tr}[X \sigma_x] $ "

JC_APO_no_roqj: Examples/JC_APO_no_roqj.cpp
	g++ Examples/JC_APO_no_roqj.cpp -o Examples/JC_APO_no_roqj.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/JC_APO_no_roqj.x
	python3 Examples/plot_JC_APO_no_roqj.py

sigma_x_sigma_z: roqj.o roqj_pop.o Examples/sigma_x_sigma_z.cpp
	g++ Examples/sigma_x_sigma_z.cpp roqj.o roqj_pop.o -o Examples/sigma_x_sigma_z.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/sigma_x_sigma_z.x
	#python3 Examples/plot_sigma_x_sigma_z.py "$$ \sigma_x - \sigma_z$$ coupling with APO, $$ \Phi^z_t(Q_x) $$, $$ \vert \Psi_0\rangle = (\vert0,0\rangle + \vert1,1\rangle)/\sqrt{2} $$ " "$$ \operatorname{tr}[X \sigma_x] $$ "

ph_cov_0_discord_2: roqj.o roqj_pop.o Examples/ph_cov_0_discord_2.cpp
	g++ Examples/ph_cov_0_discord_2.cpp roqj.o roqj_pop.o -o Examples/ph_cov_0_discord_2.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/ph_cov_0_discord_2.x
	#python3 Examples/plot_ph_cov.py "$ \Phi^0_t(Q_x) $" "$ \operatorname{tr}[X \sigma_x] $"

Pauli_channels: roqj.o roqj_pop.o Examples/Pauli_channels.cpp
	g++ Examples/Pauli_channels.cpp roqj.o roqj_pop.o -o Examples/Pauli_channels.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/Pauli_channels.x
	mv -f analytic.txt analytic_m.txt
	python3 Examples/plot_Pauli.py "Pauli channels, $$ \Phi^0_t(Q_z) $$" "$$ \operatorname{tr}[X \sigma_z] $$"