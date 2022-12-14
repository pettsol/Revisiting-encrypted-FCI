#include <cmath>
#include "acc_efci.h"

void acc_efci_encrypt(
		uint64_t *c_tr,
		uint64_t c_Px[],
		uint64_t c_P[],
		rabbit_state *state,
		const uint64_t trace,
		const uint64_t Px[],
		const uint64_t P[],
		const uint32_t dim)
{
	// Create enough pseudorandom words?
	// Number given by 1 + dim*dim + dim
	uint64_t n_bytes = sizeof(uint64_t) * (1 + dim*dim + dim);

	// Generate the bits
	uint8_t *array = new uint8_t[n_bytes];
	for (uint32_t i = 0; i < n_bytes; i++)
	{
		array[i] = 0;
	}
	rabbit_process_packet(state, array, array, n_bytes);

	// Encrypt trace
	uint64_t *word = (uint64_t*)array;

	*c_tr = trace + word[0];

	// Encrypt the information vector
	for (uint32_t i = 0; i < dim; i++)
	{
		c_Px[i] = Px[i] + word[1+i];
	}

	// Encrypt the information matrix
	for (uint32_t i = 0; i < dim*dim; i++)
	{
		c_P[i] = P[i] + word[1 + dim + i];
	}

	// Clean up
	delete[] array;
}

void acc_efci_fusion(
		uint64_t *sum_tr,
		uint64_t sum_Px[],
		uint64_t sum_P[],
		const uint64_t c_tr[],
		const uint64_t c_Px[],
		const uint64_t c_P[],
		const uint32_t dim,
		const uint32_t n_sensors)
{
	// Sum the traces
	*sum_tr = c_tr[0];
	for (uint32_t i = 1; i < n_sensors; i++)
	{
		*sum_tr += c_tr[i];
	}

	// Sum the information vector
	for (uint32_t i = 0; i < dim; i++)
	{
		sum_Px[i] = c_Px[i];
	}
	for (uint32_t i = 1; i < n_sensors; i++)
	{
		for (uint32_t j = 0; j < dim; j++)
		{
			sum_Px[j] += c_Px[i*dim + j];
		}
	}

	// Sum the information matrix
	for (uint32_t i = 0; i < dim*dim; i++)
	{
		sum_P[i] = c_P[i];
	}
	for (uint32_t i = 1; i < n_sensors; i++)
	{
		for (uint32_t j = 0; j < dim*dim; j++)
		{
			sum_P[j] += c_P[i*dim*dim + j];
		}
	}
}

void acc_efci_decrypt(
		double Px[],
		double P[],
		rabbit_state state[],
		const uint64_t sum_tr,
		const uint64_t sum_Px[],
		const uint64_t sum_P[],
		const uint64_t gamma,
		const uint32_t dim,
		const uint32_t n_sensors)
{
	// We need to generate the appropriate keystreams
	// for each sensor system
	
	// Compute the "cumulative key" in an array
	uint64_t n_bytes = sizeof(uint64_t) * (1 + dim + dim*dim);
	uint64_t *cumulative_keystream = new uint64_t[n_bytes/sizeof(uint64_t)];

	// Set array to zero
	for (uint32_t i = 0; i < n_bytes/sizeof(uint64_t); i++)
	{
		cumulative_keystream[i] = 0;
	}

	// Just iterate over all sensor systems, generate
	// the appropriate keystream and add to the cumulative
	// keystream
	uint8_t *array = new uint8_t[n_bytes];
	uint64_t *w_ptr = (uint64_t*)array;
	for (uint32_t i = 0; i < n_sensors; i++)
	{
		for (uint32_t i = 0; i < n_bytes; i++)
		{
			array[i] = 0;
		}
		rabbit_process_packet(&state[i], array, array, n_bytes);

		// Add keystream to the cumulative keystream
		for (uint32_t i = 0; i < n_bytes/sizeof(uint64_t); i++)
		{
			cumulative_keystream[i] += w_ptr[i];
		}
	}

	// Decrypt by subtracting, then decode and invert the traces
	uint64_t int_trace = sum_tr - cumulative_keystream[0];
	double traces;
	inv_rho(&traces, int_trace, gamma);
	traces = 1/traces;

	uint64_t tmp;
	// Decrypt by subtracting, then decode and multiply with the traces
	for (uint32_t i = 0; i < dim; i++)
	{
		tmp = sum_Px[i] - cumulative_keystream[1 + i];
		inv_rho(&Px[i], tmp, gamma);
		Px[i] = traces * Px[i];

	}

	for (uint32_t i = 0; i < dim*dim; i++)
	{
		tmp = sum_P[i] - cumulative_keystream[1 + dim + i];
		inv_rho(&P[i], tmp, gamma);
		P[i] = traces * P[i];
	}

	delete[] cumulative_keystream;
	delete[] array;
}

void rho(
		uint64_t *out,
		const double in,
		const uint64_t gamma)
{
	if (in < 0.0)
	{
		*out = (uint64_t)((int64_t)(in * gamma));
	}
	else
	{
		*out = (uint64_t)(in * gamma);
	}
}

void inv_rho(
		double *out,
		const uint64_t in,
		const uint64_t gamma)
{
	uint64_t half_exp = std::pow(2, 32);

	if (in < half_exp)
	{
		*out = double(in)/gamma;
	}
	else
	{
		*out = ((double)(int64_t(in)))/gamma;
	}
}
