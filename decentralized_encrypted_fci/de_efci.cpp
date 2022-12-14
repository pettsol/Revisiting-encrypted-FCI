#include <cmath>
#include "de_efci.h"

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
		const uint32_t sensor)
{
	// Encrypt trace
	ppa_encrypt(c_tr, ski, trace, *t, N, N2);
	mpz_add_ui(*t, *t, 1);

	// Encrypt the information vector
	for (uint32_t i = 0; i < dim; i++)
	{
		ppa_encrypt(c_Px[i][sensor], ski, Px[i][sensor], *t, N, N2);
		mpz_add_ui(*t, *t, 1);
	}

	// Encrypt the information matrix
	for (uint32_t i = 0; i < dim*dim; i++)
	{
		ppa_encrypt(c_P[i][sensor], ski, P[i][sensor], *t, N, N2);
		mpz_add_ui(*t, *t, 1);
	}
}

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
		const uint32_t n_sensors)
{
	// Sum and decrypt the traces
	mpz_t aggregated_trace;
	mpz_init(aggregated_trace);
	ppa_aggregate_decrypt(aggregated_trace, sk0, *t, c_tr, N, N2, n_sensors);
	mpz_add_ui(*t, *t, 1);

	// Sum the information vector
	mpz_t fused_Px[dim];
	for (uint32_t i = 0; i < dim; i++)
	{
		mpz_init(fused_Px[i]);
		ppa_aggregate_decrypt(fused_Px[i], sk0, *t, c_Px[i], N, N2, n_sensors);
		mpz_add_ui(*t, *t, 1);
	}

	// Sum the information matrix
	mpz_t fused_P[dim*dim];
	for (uint32_t i = 0; i < dim*dim; i++)
	{
		mpz_init(fused_P[i]);
		ppa_aggregate_decrypt(fused_P[i], sk0, *t, c_P[i], N, N2, n_sensors);
		mpz_add_ui(*t, *t, 1);
	}

	// Recover corresponding doubles
	double trace;
	inv_rho(&trace, aggregated_trace, gamma, N);
	trace = 1/trace;

	for (uint32_t i = 0; i < dim; i++)
	{
		inv_rho(&Px[i], fused_Px[i], gamma, N);
		Px[i] = trace*Px[i];
	}

	for (uint32_t i = 0; i < dim*dim; i++)
	{
		inv_rho(&P[i], fused_P[i], gamma, N);
		P[i] = trace*P[i];
	}
}

void rho(
		mpz_t out,
		const mpf_t in,
		const mpz_t gamma,
		const mpz_t ptspace)
{
        mpf_t tmp;
        mpf_init(tmp);
        mpf_set_z(tmp, gamma);
        mpf_mul(tmp, in, tmp);

        mpz_set_f(out, tmp);
        mpz_mod(out, out, ptspace);
}

void inv_rho(
		double *out,
		const mpz_t in,
		const mpz_t gamma,
		const mpz_t ptspace)
{
        mpz_t halfsize;
        mpz_init(halfsize);

        mpz_div_ui(halfsize, ptspace, 2);
        mpz_t test;
        mpz_init(test);

        mpz_sub(test, in, halfsize);

        mpf_t gamma_f;
        mpf_init(gamma_f);

        mpf_set_z(gamma_f, gamma);

        if (mpz_sgn(test) != -1)
        {
                // negative number
                mpz_sub(test, in, ptspace);
        } else {
                mpz_set(test, in);
        }
	mpf_t tmp_out;
	mpf_init(tmp_out);
        mpf_set_z(tmp_out, test);
        mpf_div(tmp_out, tmp_out, gamma_f);
	*out = mpf_get_d(tmp_out);
}
