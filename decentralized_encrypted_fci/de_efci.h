#ifndef DE_EFCI_H
#define DE_EFCI_H

#include "privacy_preserving_aggregation/ppa.h"

void de_efci_encrypt(
		mpz_t c_tr,
		mpz_t **c_Px,
		mpz_t **c_P,
		mpz_t *t,
		const mpz_t trace,
		mpz_t **Px,
		mpz_t **P,
		const mpz_t ski,
		const mpz_t N,
		const mpz_t N2,
		const uint32_t dim,
		const uint32_t sensor);

void de_efci_fuse_and_decrypt(
		double Px[],
		double P[],
		mpz_t *t,
		const mpz_t c_tr[],
		mpz_t **c_Px,
		mpz_t **c_P,
		const mpz_t sk0,
		const mpz_t N,
		const mpz_t N2,
		const mpz_t gamma,
		const uint32_t dim,
		const uint32_t n_sensors);

void rho(
		mpz_t out,
		const mpf_t in,
		const mpz_t gamma,
		const mpz_t ptspace);

void inv_rho(
		double *out,
		const mpz_t in,
		const mpz_t gamma,
		const mpz_t ptspace);
#endif
