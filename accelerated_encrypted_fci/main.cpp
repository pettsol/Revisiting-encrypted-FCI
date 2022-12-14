#include <fstream>
#include <iostream>
#include <unistd.h>
#include <chrono>
#include <cmath>

#include "acc_efci.h"

int main()
{
	uint32_t timesteps = 150;
	uint32_t n_datasets = 1;
	uint32_t KEYSIZE = 16;
	uint32_t keysize;

	uint32_t dim = 4;
	//uint32_t n_sensors = 2;

	uint32_t sensor_array[] = {2, 4, 6, 8, 10, 15};
	size_t total(sizeof(sensor_array) / sizeof(uint32_t));

        uint32_t keysize_array[] = {128};
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

        sprintf(key_buf, "latency/keysize_%u/decryption_latency.csv",
                        keysize_array[n_keysize]);
        std::string decryption_latency = key_buf;

        std::ofstream Encryption_latency_file, Fusion_latency_file, Decryption_latency_file;

        Encryption_latency_file.open(encrypt_latency);
        Fusion_latency_file.open(fusion_latency);
        Decryption_latency_file.open(decryption_latency);


	// WE ITERATE OVER ALL NUMBERS OF SENSORS
	for (uint32_t number = 0; number < total; number++)
	{
	
	uint32_t n_sensors = sensor_array[number];

	// Choose keys and initialize N sensor systems?
	uint8_t *keys = new uint8_t[n_sensors * KEYSIZE];
	uint8_t *iv = new uint8_t[8];

	// Pick keys that are not the same to verify correct
	// implementation. Initialize the states
	uint64_t *w_ptr = (uint64_t*)keys;
	// We need two sets of states. One for encryption
	// and one for decryption.
	rabbit_state *enc_states = new rabbit_state[n_sensors];
	rabbit_state *dec_states = new rabbit_state[n_sensors];
	for (uint64_t i = 0; i < n_sensors; i++)
	{
		w_ptr[2*i] += i;
		rabbit_load_key(&enc_states[i], &keys[i * KEYSIZE]);
		rabbit_load_key(&dec_states[i], &keys[i * KEYSIZE]);
		rabbit_load_iv(&enc_states[i], iv);
		rabbit_load_iv(&dec_states[i], iv);
	}

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

        uint32_t count;
        count = dim*dim*n_sensors*timesteps;

        double *Pm_double_array = new double[count];
        double *P_double_array = new double[count];

        for (uint32_t i = 0; i < count; i++)
        {
                Pm_file >> Pm_double_array[i];
                P_file >> P_double_array[i];
        }

        count = count / dim;
        double *Px_double_array = new double[count];
        for (uint32_t i = 0; i < count; i++)
        {
                Px_file >> Px_double_array[i];
        }

	double *traces = new double[n_sensors*timesteps];

	// Map the arrays to unsigned integers
        uint64_t gamma = std::pow(2, 40);

	uint64_t *Trace = new uint64_t[n_sensors*timesteps];
        uint64_t *P = new uint64_t[count*dim];
        uint64_t *Px = new uint64_t[count];

	double *P0 = new double[dim*dim];
	double *P0x0 = new double[dim];

	double tmp;
	// Compute the inverse of the traces for each sensor for each timestep
	for (uint32_t i = 0; i < timesteps; i++)
	{
		for (uint32_t j = 0; j < n_sensors; j++)
		{
			traces[i*n_sensors + j] = 0;
			for (uint32_t k = 0; k < dim; k++)
			{
				traces[i*n_sensors + j] += Pm_double_array[
					i*dim*dim*n_sensors + 
					j*dim*dim + k*dim + k];
			}
			traces[i*n_sensors + j] = 1/traces[i*n_sensors + j];
			// Rho
			rho(&Trace[i*n_sensors + j], traces[i*n_sensors + j],
				       	gamma);
		}
	}

	// Multiply all matrices with the inverse trace and then map to plaintexts
	for (uint32_t i = 0; i < timesteps; i++)
	{
		for (uint32_t j = 0; j < n_sensors; j++)
		{
			for (uint32_t k = 0; k < dim*dim; k++)
			{
				tmp = P_double_array[i*dim*dim*n_sensors +
				       	j*dim*dim + k] *
				       traces[i*n_sensors + j];
				rho(&P[i*dim*dim*n_sensors + j*dim*dim + k], tmp,
					       	gamma);
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
				tmp = Px_double_array[i*dim*n_sensors + j*dim + k] *
				       traces[i*n_sensors + j];
				rho(&Px[i*dim*n_sensors + j*dim + k], tmp,
					       	gamma);
			}
		}
	}

	uint64_t *c_tr = new uint64_t[n_sensors];
	uint64_t *c_Px = new uint64_t[dim*n_sensors];
	uint64_t *c_P = new uint64_t[dim*dim*n_sensors];
	uint64_t *sum_P = new uint64_t[dim*dim];
	uint64_t *sum_Px = new uint64_t[dim];
	uint64_t sum_tr = 0;

	for (uint32_t i = 0; i < n_sensors; i++)
	{
		c_tr[i] = 0;
		for (uint32_t j = 0; j < dim; j++)
		{
			c_Px[dim*i + j] = 0;
			for (uint32_t k = 0; k < dim; k++)
			{
				c_P[i*dim*dim + j*dim + k] = 0;
			}
		}	
	}

	for (uint32_t i = 0; i < dim; i++)
	{
		sum_Px[i] = 0;
		for (uint32_t j = 0; j < dim; j++)
		{
			sum_P[i*dim + j] = 0;
		}
	}

	// Iterate over timesteps
	for (uint32_t t = 0; t < timesteps; t++)
	{
		// Encrypt the sensor data
		auto start = std::chrono::high_resolution_clock::now();
		for (uint32_t i = 0; i < n_sensors; i++)
		{
                        // Start clock
                        start = std::chrono::high_resolution_clock::now();
			acc_efci_encrypt(&c_tr[i], &c_Px[i*dim], &c_P[i*dim*dim],
					&enc_states[i], Trace[t*n_sensors + i], 
					&Px[t*dim*n_sensors + i*dim],
					&P[t*n_sensors*dim*dim + i*dim*dim],
					dim);
		}
                // Only write the most recent encryption latency. Writing all is not required.
                // Stop clock
                auto stop = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>
                                (stop - start);

                Encryption_latency_file << duration.count()  << " ";


		// Proceed by fusing the encrypted data
                std::cout << "*** FUSING ENCRYPTED DATA FROM " << sensor_array[number]  << " SENSORS AT TIMESTEP " << t+1 << " IN DATASET " << dataset+1 << " ***" << std::endl;
                // Start clock
                start = std::chrono::high_resolution_clock::now();
		// Fuse the data
		acc_efci_fusion(&sum_tr, sum_Px, sum_P, c_tr, 
				c_Px, c_P, dim, n_sensors);
                // Stop clock
                stop = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
                Fusion_latency_file << duration.count() << " ";


		// Start clock
                start = std::chrono::high_resolution_clock::now();
		// Decrypt the data
		acc_efci_decrypt(P0x0, P0, dec_states, sum_tr, sum_Px, sum_P,
				gamma, dim, n_sensors);

                // Stop clock
                stop = std::chrono::high_resolution_clock::now();

                duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
                Decryption_latency_file << duration.count() << " ";


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
       	delete[] Pm_double_array;
        delete[] P_double_array;
        delete[] Px_double_array;
        delete[] P;
        delete[] Px;
	delete[] Trace;
        delete[] c_tr;
	delete[] c_P;
	delete[] c_Px;
	delete[] traces;
        delete[] sum_Px;
        delete[] sum_P;
        delete[] P0;
        delete[] P0x0;
	}
	// Add a newline
        Encryption_latency_file << std::endl;
        Fusion_latency_file << std::endl;
        Decryption_latency_file << std::endl;
        
	// Clean up
	delete[] keys;
	delete[] iv;
	delete[] enc_states;
	delete[] dec_states;
	}
        // Close the files
        Encryption_latency_file.close();
        Fusion_latency_file.close();
        Decryption_latency_file.close();
	}

	std::cout << "FUSION OF ALL DATASETS COMPLETE" << std::endl;
}
