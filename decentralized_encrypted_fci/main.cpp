#include <fstream>
#include <iostream>
#include <unistd.h>
#include <chrono>
#include <cmath>

#include "de_efci.h"

int main()
{
	const uint32_t timesteps = 150;
	const uint32_t n_datasets = 100;
	uint32_t keysize;

	const uint32_t dim = 4;
	//uint32_t n_sensors = 2;

	uint32_t sensor_array[] = {2, 4, 6, 8, 10, 15};
	size_t total(sizeof(sensor_array) / sizeof(uint32_t));

        uint32_t keysize_array[] = {1024};
        size_t total_keysizes(sizeof(keysize_array) / sizeof(uint32_t));

        // We iterate over all keysizes
        for (uint32_t n_keysize = 0; n_keysize < total_keysizes; n_keysize++)
        {

        keysize = keysize_array[n_keysize];
        char key_buf[100];
        sprintf(key_buf, "latency/keysize_%u/sensor_encryption_latency.csv",
                        keysize_array[n_keysize]);
        std::string encrypt_latency = key_buf;

        sprintf(key_buf, "latency/keysize_%u/encrypted_fusion_latency.csv",
                        keysize_array[n_keysize]);
        std::string fusion_latency = key_buf;

        std::ofstream Encryption_latency_file, Fusion_latency_file;

        Encryption_latency_file.open(encrypt_latency);
        Fusion_latency_file.open(fusion_latency);


	// WE ITERATE OVER ALL NUMBERS OF SENSORS
	for (uint32_t number = 0; number < total; number++)
	{
	
	const uint32_t n_sensors = sensor_array[number];

	// Initialize keys for N sensor systems
	mpz_t *keys = new mpz_t[n_sensors+1];

	for (uint64_t i = 0; i < n_sensors+1; i++)
	{
		mpz_init(keys[i]);
	}

	// Generate the keys and the factoring modulus
	mpz_t N, N2;
	mpz_init(N);
	mpz_init(N2);
	ppa_setup(N, keys, keysize, n_sensors);
	mpz_mul(N2, N, N);

	// WE ITERATE THROUGH ALL DATASETS

	for (uint32_t dataset = 0; dataset < n_datasets; dataset++)
	{
	// We need to extract information from the files and store in 3 arrays
        char buf[100];
        sprintf(buf, "datasets/n_sensors=%u/covariance_%u.csv", sensor_array[number], dataset);
        std::string covariance_string = buf;
        sprintf(buf, "datasets/n_sensors=%u/inv_covariance_%u.csv", sensor_array[number], dataset);
        std::string inv_covariance_string = buf;
        sprintf(buf, "datasets/n_sensors=%u/inv_covariance_x_mean_%u.csv", sensor_array[number], dataset);
        std::string inv_covariance_x_mean_string = buf;

        std::ifstream Pm_file, P_file, Px_file;
        Pm_file.open(covariance_string);
        P_file.open(inv_covariance_string);
        Px_file.open(inv_covariance_x_mean_string);

        // We need to write the fused P0 and P0x0 to files
        sprintf(buf, "output/n_sensors=%u/inv_covariance_%u.csv", sensor_array[number], dataset);
        std::string inv_covariance_string_out = buf;
        sprintf(buf, "output/n_sensors=%u/inv_covariance_x_mean_%u.csv", sensor_array[number], dataset);
        std::string inv_covariance_x_mean_string_out = buf;

        std::ofstream P0_file(inv_covariance_string_out);
        std::ofstream P0x0_file(inv_covariance_x_mean_string_out);


	// EDIT TO EXTRACT TO NEW DATASTRUCTURE
        const uint32_t count = dim*dim*timesteps;

        double **Pm_double_array = new double*[count];
        double **P_double_array = new double*[count];
	for (uint32_t i = 0; i < count; i++)
	{
		Pm_double_array[i] = new double[n_sensors];
		P_double_array[i] = new double[n_sensors];
	}

	const uint32_t count_v = count / dim;
        double **Px_double_array = new double*[count_v];

	for (uint32_t i = 0; i < count_v; i++)
	{
		Px_double_array[i] = new double[n_sensors];
	}

	// Iterate over all timesteps, over all sensors, over all elements
	for (uint32_t i = 0; i < timesteps; i++)
	{
		for (uint32_t j = 0; j < n_sensors; j++)
		{
			// Extract Pm and P
			for (uint32_t k = 0; k < dim*dim; k++)
			{
				Pm_file >> Pm_double_array[i*dim*dim + k][j];
				P_file >> P_double_array[i*dim*dim + k][j];
			}

			// Extract Px
			for (uint32_t k = 0; k < dim; k++)
			{
				Px_file >> Px_double_array[i*dim + k][j];
			}
		}
	}

	Pm_file.close();
	P_file.close();
	Px_file.close();

	double *traces = new double[n_sensors*timesteps];

	// Map the arrays to unsigned integers
        mpz_t gamma;// = std::pow(2, 40);
	mpz_init(gamma);
	mpz_ui_pow_ui(gamma, 2, 40);

	mpz_t *Trace = new mpz_t[n_sensors*timesteps];
        mpz_t **P = new mpz_t*[count];
        mpz_t **Px = new mpz_t*[count_v];

	for (uint32_t i = 0; i < count; i++)
	{
		P[i] = new mpz_t[n_sensors];
	}

	for (uint32_t i = 0; i < count_v; i++)
	{
		Px[i] = new mpz_t[n_sensors];
	}

	double *P0 = new double[dim*dim];
	double *P0x0 = new double[dim];

	mpf_t input;
	mpf_init(input);

	double tmp;
	// Compute the inverse of the traces for each sensor for each timestep
	for (uint32_t i = 0; i < timesteps; i++)
	{
		for (uint32_t j = 0; j < n_sensors; j++)
		{
			mpz_init_set_ui(Trace[i*n_sensors + j], 0);
			traces[i*n_sensors + j] = 0;
			for (uint32_t k = 0; k < dim; k++)
			{
				traces[i*n_sensors + j] += Pm_double_array[i*dim*dim + k*dim + k][j];
			}
			traces[i*n_sensors + j] = 1/traces[i*n_sensors + j];
			// Map to valid plaintext
			mpf_set_d(input, traces[i*n_sensors + j]);
			rho(Trace[i*n_sensors + j], input,
				       	gamma, N);
		}
	}

	// Multiply all matrices with the inverse trace and then map to plaintexts
	for (uint32_t i = 0; i < timesteps; i++)
	{
		for (uint32_t j = 0; j < n_sensors; j++)
		{
			for (uint32_t k = 0; k < dim*dim; k++)
			{
				mpz_init(P[i*dim*dim + k][j]);
				tmp = P_double_array[i*dim*dim + k][j] *
				       traces[i*n_sensors + j];
				// Map to valid plaintext
				mpf_set_d(input, tmp);
				rho(P[i*dim*dim + k][j], input,
					  gamma, N);
			}
		}
	}

	// Multiply all vectors with the inverse trace and then map to plaintexts
	for (uint32_t i = 0; i < timesteps; i++)
	{
		for (uint32_t j = 0; j < n_sensors; j++)
		{
			for (uint32_t k = 0; k < dim; k++)
			{
				mpz_init(Px[i*dim + k][j]);
				tmp = Px_double_array[i*dim + k][j] *
				       traces[i*n_sensors + j];
				// Map to plaintext
				mpf_set_d(input, tmp);
				rho(Px[i*dim + k][j], input,
					gamma, N);
			}
		}
	}

	mpz_t *c_tr = new mpz_t[n_sensors];
	mpz_t **c_Px = new mpz_t*[dim];
	mpz_t **c_P = new mpz_t*[dim*dim];

	for (uint32_t i = 0; i < dim*dim; i++)
	{
		c_P[i] = new mpz_t[n_sensors];
	}

	for (uint32_t i = 0; i < dim; i++)
	{
		c_Px[i] = new mpz_t[n_sensors];
	}

	for (uint32_t i = 0; i < n_sensors; i++)
	{
		mpz_init_set_ui(c_tr[i], 0);
		for (uint32_t j = 0; j < dim; j++)
		{
			mpz_init_set_ui(c_Px[j][i], 0);
			for (uint32_t k = 0; k < dim; k++)
			{
				mpz_init_set_ui(c_P[j*dim + k][i], 0);
			}
		}	
	}

	mpz_t timestep;
	mpz_init_set_ui(timestep, 0);
	// Iterate over timesteps
	for (uint32_t t = 0; t < timesteps; t++)
	{
		// Encrypt the sensor data
		auto start = std::chrono::high_resolution_clock::now();
		for (uint32_t i = 0; i < n_sensors; i++)
		{
			mpz_t old_timestep;
		        mpz_init_set(old_timestep, timestep);
                        // Start clock
                        start = std::chrono::high_resolution_clock::now();
			de_efci_encrypt(c_tr[i], c_Px, c_P, &timestep,
					Trace[t*n_sensors + i], &Px[t*dim], &P[t*dim*dim], keys[i+1], N, N2, dim, i);
			mpz_set(timestep, old_timestep);
		}
                // Only write the most recent encryption latency. Writing all is not required.
                // Stop clock
                auto stop = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>
                                (stop - start);

                Encryption_latency_file << duration.count()  << " ";


		// Proceed by fusing and decrypting the encrypted data
                std::cout << "*** FUSING ENCRYPTED DATA FROM " << sensor_array[number]  << " SENSORS AT TIMESTEP " << t+1 << " IN DATASET " << dataset+1 << " ***" << std::endl;
                // Start clock
                start = std::chrono::high_resolution_clock::now();
		// Fuse and decrypt the data
		de_efci_fuse_and_decrypt(P0x0, P0, &timestep, c_tr, 
				c_Px, c_P, keys[0], N, N2, gamma, dim, n_sensors);
                // Stop clock
                stop = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
                Fusion_latency_file << duration.count() << " ";


		// Write to file
		for (uint32_t i = 0; i < dim*dim; i++)
		{
			P0_file << P0[i];
			P0_file << " ";
		}
		
		// Timesteps separated by newline
		P0_file << std::endl;

		for (uint32_t i = 0; i < dim; i++)
		{
			P0x0_file << P0x0[i];
			P0x0_file << " ";
		}
		P0x0_file << std::endl;
	}

	// Close the files
	P0_file.close();
	P0x0_file.close();

	// Clean up
	for (uint32_t i = 0; i < count; i++)
	{
		delete[] Pm_double_array[i];
		delete[] P_double_array[i];
		delete[] P[i];
	}
	delete[] Pm_double_array;
	delete[] P_double_array;
	delete[] P;

	for (uint32_t i = 0; i < count_v; i++)
	{
		delete[] Px_double_array[i];
		delete[] Px[i];
	}
	delete[] Px_double_array;
	delete[] Px;

	for (uint32_t i = 0; i < dim*dim; i++)
	{
		delete[] c_P[i];
	}
	delete[] c_P;

	for (uint32_t i = 0; i < dim; i++)
	{
		delete[] c_Px[i];
	}
	delete[] c_Px;

	delete[] Trace;
        delete[] c_tr;
	delete[] traces;
        delete[] P0;
        delete[] P0x0;
	}
	// Add a newline
        Encryption_latency_file << std::endl;
        Fusion_latency_file << std::endl;
        
	// Clean up
	delete[] keys;
	}
        // Close the files
        Encryption_latency_file.close();
        Fusion_latency_file.close();
	}

	std::cout << "FUSION OF ALL DATASETS COMPLETE" << std::endl;
}
