
#include "Rabbit/rabbit.h"

void acc_efci_encrypt(
		uint64_t *c_tr,
		uint64_t c_Px[],
		uint64_t c_P[],
		rabbit_state *state,
		const uint64_t trace,
		const uint64_t Px[],
		const uint64_t P[],
		const uint32_t dim);

void acc_efci_fusion(
		uint64_t *sum_tr,
		uint64_t sum_Px[],
		uint64_t sum_P[],
		const uint64_t c_tr[],
		const uint64_t c_Px[],
		const uint64_t c_P[],
		const uint32_t dim,
		const uint32_t n_sensors);

void acc_efci_decrypt(
		double Px[],
		double P[],
		rabbit_state state[],
		const uint64_t sum_tr,
		const uint64_t sum_Px[],
		const uint64_t sum_P[],
		const uint64_t gamma,
		const uint32_t dim,
		const uint32_t n_sensors);

void rho(
		uint64_t *out,
		const double in,
		const uint64_t gamma);

void inv_rho(
		double *out,
		const uint64_t in,
		const uint64_t gamma);
