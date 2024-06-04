#ifndef __HGSW_FUNCS__
#define __HGSW_FUNCS__

#include "lattice/dgsampling.h"
#include "rng/sampler.h"

using namespace std;
using namespace lbcrypto;
using namespace hgsw_params;

namespace hgsw_funcs {
    static inline uint64_t rp(__uint128_t a, __uint128_t b, const uint64_t q, const __uint128_t r) // Plantard's modular reduction mod p
    {
	    __uint128_t c;
	    c = a * b * r;
	    return (((c >> prm::rp_shift) + 1) * q) >> prm::rp_shift;
    }

    static inline uint64_t rp_const(__uint128_t a, __uint128_t b, const uint64_t q) // auxiliary function to perform NTT on plantard form
    {
	    __uint128_t c;
	    c = a * b;
	    return (((c >> prm::rp_shift) + 1) * q) >> prm::rp_shift;
    }

    static inline uint64_t con_sub(const uint64_t x, const uint64_t q) // x - q if x >= q
    {
	    return x - ((-(1 ^ ((x - q) >> 63))) & q);
    }

    static inline uint64_t con_add(const uint64_t x, const uint64_t q) // x + q if x <= 0
    {
	    return x + ((-(x >> 63)) & q);
    }

    static inline uint16_t rp(uint32_t a) // auxiliary function used in multiplication mod p in Plantard form  
    {
	    uint32_t c;
	    c = a * 226050911;
	    return (((c >> 16) + 1) * prm::p) >> 16;
    }

    static void rp_mulp(uint16_t *out, const uint16_t *a, const uint16_t *b) // multiplication mod p in Plantard form
    {
	    size_t rp_len = prm::d / 2;
	    int16_t c[prm::d * 2] = {0}, z[prm::d] = {0};
	    int16_t a1[rp_len], b1[rp_len];
	    uint16_t i, j;
	    int16_t x;
	    const uint16_t p0 = prm::p * 546;

	    /* Layer 1 Karatsuba+schoolbook multiplication */
	    for (i = 0; i < rp_len; i++) {
		    a1[i] = a[i] - a[i + rp_len];
		    b1[i] = b[i + rp_len] - b[i];
	    }
	    for (i = 0; i < rp_len; i++) {
		    for (j = 0; j < rp_len; j++) {
			    c[i + j] = c[i + j] + a[i] * b[j];
			    c[i + j + prm::d] = c[i + j + prm::d] + a[i + rp_len] * b[j + rp_len];
			    z[i + j] = z[i + j] + a1[i] * b1[j];
		    }
	    }
	    for (i = 0; i < rp_len; i++) {
		    x = c[i + rp_len] + c[i + prm::d];
		    c[i + rp_len] = x + c[i] + z[i];
		    c[i + prm::d] = x + c[i + rp_len * 3] + z[i + rp_len];
	    }

	    /* mod p mod x^d+1 */
	    for (i = 0; i < prm::d; i++) {
		    out[i] = rp(c[i] - c[i + prm::d] + p0);
	    }
    }

    /* Cooley-Tukey NTT */
    void NTT(uint64_t *a, const uint64_t q, const __uint128_t conv_fact, const __uint128_t *root) // Number Theoretic Transform working on Plantard form
    {
	    uint64_t t = prm::d >> 1;
	    uint64_t m, i, j;
	    uint64_t j1, j2;
	    uint64_t x, u;
	    
	    for (j = 0; j < t; j++) {
		    x = rp_const(a[j], conv_fact, q);
		    u = rp_const(a[j + t], root[1], q);
		    a[j] = con_sub(x + u, q);
		    a[j + t] = con_add(x - u, q);
	    }
	    
	    for (m = 2; m < (prm::d >> 1); m <<= 1) {
		    t >>= 1;
		    for (i = 0; i < m; i++) {
			    j1 = (i * t) << 1;
			    j2 = j1 + t - 1;
			    
			    for (j = j1; j <= j2; j++) {
				    x = a[j];
				    u = rp_const(a[j + t], root[m + i], q);
				    a[j] = con_sub(x + u, q);
				    a[j + t] = con_add(x - u, q);
			    }
		    }
	    }

	    for (i = 0; i < m; i++) {
		    j = i << 1;
		    
		    x = a[j];
		    u = rp_const(a[j + 1], root[m + i], q);
		    a[j] = con_sub(x + u, q);
		    a[j + 1] = con_add(x - u, q);
	    }
    }

    /* Gentleman-Sande INTT */
    void INTT(uint64_t *a, const uint64_t q, const __uint128_t inv_n, const __uint128_t *root_inv) // Inverse Number Theoretic Transform working on Plantard form
    {
	    uint64_t t = 2;
	    uint64_t m, i, j;
	    uint64_t j1, j2;
	    uint64_t h = prm::d >> 1;
	    uint64_t x, y;
	    
	    for (i = 0; i < h; i++) {
		    j = i << 1;
		    
		    x = con_sub(a[j] + a[j + 1], q);
		    y = con_add(a[j] - a[j + 1], q);
		    a[j] = x;
		    a[j + 1] = rp_const(y, root_inv[h + i], q);
	    }
	    
	    for (m = prm::d >> 1; m > 2; m >>= 1) {
		    j1 = 0;
		    h >>= 1;

		    for (i = 0; i < h; i++) {
			    j2 = j1 + t - 1;
			    
			    for (j = j1; j <= j2; j++) {
				    x = con_sub(a[j] + a[j + t], q);
				    y = con_add(a[j] - a[j + t], q);
				    a[j] = x;
				    a[j + t] = rp_const(y, root_inv[h + i], q);
			    }
			    j1 += t << 1;
		    }
		    t <<= 1;
	    }

	    for (j = 0; j < t; j++) {
		    x = con_sub(a[j] + a[j + t], q);
		    y = con_add(a[j] - a[j + t], q);
		    a[j] = rp_const(x, inv_n, q);
		    a[j + t] = rp_const(y, root_inv[1], q);
	    }
    }

    template <typename T, typename prm, uint32_t ROW, uint32_t COLUMN> class mat_hgsw {
    public:
        
        vector<array<T, COLUMN>> mat;

        mat_hgsw() : mat(ROW) {}

        mat_hgsw(const mat_hgsw &other) = default;

        inline array<T, COLUMN> &operator[](const uint64_t &i) {
            return this->mat[i];
        }

        inline const array<T, COLUMN> &
        operator[](const uint64_t &i) const {
            return this->mat[i];
        }

        inline mat_hgsw &gen_gT_p() { // generate gT vector in NTT form
            for (uint32_t j = 0; j < COLUMN; j++) {
                this->mat[0][j][0] = (uint64_t)pow(prm::beta, j) % prm::p;
                NTT(this->mat[0][j], prm::qbar, prm::rp_conv_fact_qbar, prm::root_qbar);
            }
            return *this;
        }

        inline mat_hgsw &gen_qi() { // generate an uniformly distributed random qi vector
            uint8_t seed[prm::d];
	        aes256ctr_ctx ctx;
	        randombytes(seed, prm::d);
	        aes256ctr_init(&ctx, seed, 0);
            for (uint32_t j = 0; j < COLUMN; j++) {
                uni_sampler_zzp_1d(this->mat[0][j], prm::f, &ctx);
            }
            return *this;
        }

        inline mat_hgsw &gen_a() { // generate scaling factor a with uniform distribution
            uint8_t seed[prm::d];
	        aes256ctr_ctx ctx;
	        randombytes(seed, prm::d);
	        aes256ctr_init(&ctx, seed, 0);
            uni_sampler_zzp_1d(this->mat[0][0], prm::f, &ctx);
            return *this;
        }
        
        inline mat_hgsw &get_f_arr() { // fill f_arr which is used in the polynomial modular reduction in CRT Encode
            for (uint32_t i = 0; i < COLUMN; i++) {
                for (uint32_t j = 0; j < prm::d; j++) {
                    this->mat[0][i][j] = 0;
                }
            }    
            // 1 = a * f1 + b * f2 
            this->mat[0][0][16] = 1; this->mat[0][0][8] = 6;  this->mat[0][0][0] = 18; // f1 
            this->mat[0][1][16] = 1; this->mat[0][1][8] = 13; this->mat[0][1][0] = 18; // f2
            this->mat[0][2][8] = 8;  this->mat[0][2][0] = 9; // a
            this->mat[0][3][8] = 11; this->mat[0][3][0] = 9; // b
            return *this;
        }

        inline mat_hgsw &print_matrix() { // function to print a matrin in the form of mat_hgsw
            for (uint32_t i = 0; i < ROW; i++) {
                for (uint32_t j = 0; j < COLUMN; j++) {
                    for (uint32_t k = 0; k < prm::d; k++) {
                        cout << this->mat[i][j][k] << " ";   
                    }
                    printf("\n");
                }
                printf("\n\n");
            }        
            return *this;
        }

        inline mat_hgsw &sample_matrix_uniform_p() { // output a uniformly distributed random matrix mod p 
            uint8_t seed[32];
            aes256ctr_ctx ctx;           
            randombytes(seed, prm::d);
            aes256ctr_init(&ctx, seed, 0);              
            for (uint32_t i = 0; i < ROW; i++) {
                for (uint32_t j = 0; j < COLUMN; j++) {
                    uni_sampler_zzp_1d(this->mat[i][j], prm::d, &ctx);
                }
            }
            return *this;
        }

        inline mat_hgsw &sample_matrix_gaussian(DiscreteGaussianGenerator dgg) {  // output a matrix mod q generated via discrete gaussian distribution
            for (uint32_t i = 0; i < ROW; i++) {
                for (uint32_t j = 0; j < COLUMN; j++) {
                    for (uint32_t k = 0; k < prm::d; k++) {
                        this->mat[i][j][k] = (dgg.GenerateIntegerKarney(0, prm::sig_s) + prm::q) % prm::q;
                    }
                }
            }        
            return *this;
        }

        template <uint32_t C2>
        inline mat_hgsw<T, prm, ROW, C2>
        operator*(const mat_hgsw<T, prm, COLUMN, C2> &o) const { // multiply given two matrices in NTT form mod qbar
            mat_hgsw<T, prm, ROW, C2> res;
            for (uint32_t i = 0; i < ROW; i++) {
                for (uint32_t j = 0; j < C2; j++) {
                    for (uint32_t k = 0; k < COLUMN; k++) {
                        T c;
                        for (uint32_t l = 0; l < prm::d; l++) {
                            c[l] = rp(this->mat[i][k][l], o[k][j][l], prm::qbar, prm::rp_r_qbar);
                            res[i][j][l] = (res[i][j][l] + c[l]) % prm::qbar;
                        }
                    }
                }
            }
            return res;
        }
        
        inline mat_hgsw &operator+=(const mat_hgsw &o) { // add given two matrices in NTT form mod qbar
            for (uint32_t i = 0; i < ROW; i++) {
                for (uint32_t j = 0; j < COLUMN; j++) {
                    for (uint32_t k = 0; k < prm::d; k++) {
                        this->mat[i][j][k] = (this->mat[i][j][k] + o[i][j][k]) % prm::qbar;
                    }
                }
            }
            return *this;
        }

        inline mat_hgsw &add_mod(const mat_hgsw &o, uint64_t modulus) { // add given two matrices according to given modulo
            for (uint32_t i = 0; i < ROW; i++) {
                for (uint32_t j = 0; j < COLUMN; j++) {
                    for (uint32_t k = 0; k < prm::d; k++) {
                        this->mat[i][j][k] = (this->mat[i][j][k] + o[i][j][k]) % modulus;
                    }
                }
            }
            return *this;
        }
        
        template <uint32_t C2>
        inline mat_hgsw &mul_mod_qprop(const mat_hgsw<T, prm, ROW, C2> &mat1, const mat_hgsw<T, prm, C2, COLUMN> &mat2) { // multiply given two matrices in NTT form mod qpr/p 
            for (uint32_t i = 0; i < ROW; i++) {
                for (uint32_t j = 0; j < COLUMN; j++) {
                    for (uint32_t k = 0; k < C2; k++) {
                        T c;
                        for (uint32_t l = 0; l < prm::d; l++) {
                            c[l] = rp(mat1[i][k][l], mat2[k][j][l], prm::qprop, prm::rp_r_qprop);
                            this->mat[i][j][l] = (this->mat[i][j][l] + c[l]) % prm::qprop;
                        }
                    }
                }
            } 
            return *this;       
        } 
        
        template <uint32_t C2>
        inline mat_hgsw &mul_modp(const mat_hgsw<T, prm, ROW, C2> &mat1, const mat_hgsw<T, prm, C2, COLUMN> &mat2) { // multiply given two matrices mod p
            for (uint32_t i = 0; i < ROW; i++) {
                for (uint32_t j = 0; j < COLUMN; j++) {
                    for (uint32_t k = 0; k < C2; k++) {
                        T c;
                        rp_mulp(c, mat1[i][k], mat2[k][j]);
                        for (uint32_t l = 0; l < prm::d; l++) {
                           this->mat[i][j][l] = (this->mat[i][j][l] + c[l]) % prm::p;
                        }
                    }
                }
            } 
            return *this;       
        }  

        inline mat_hgsw &compute_Hit_p(const mat_hgsw<T, prm, 1, COLUMN> &mat1, const mat_hgsw<T, prm, 1, ROW> &mat2) { // special function to perform the operation at Appx.F Alg.2 Line 6-7 
            for (uint32_t i = 0; i < COLUMN; i++) { 
                for (uint32_t j = 0; j < ROW; j++) { 
                    for (uint32_t k = 0; k < prm::d; k++) {
                        this->mat[j][i][k] = rp(mat1[0][i][k], mat2[0][j][0], prm::qbar, prm::rp_r_qbar);
                        this->mat[j][i][k] = rp(prm::qbmp, this->mat[j][i][k], prm::qbar, prm::rp_r_qbar);
                    }
                }
            } 
            return *this;       
        }      
        
        inline mat_hgsw &scale_convert(const mat_hgsw<T, prm, ROW, COLUMN> &H_bar) { // special function to perform the operation at Appx.F Alg.4 Line 6
            for (uint32_t i = 0; i < COLUMN; i++) {
                for (uint32_t j = 0; j < prm::d; j++){ 
                    double dblval = H_bar[0][i][j];
                    this->mat[0][i][j] = (uint16_t)round(prm::poqpr * dblval) % prm::p;
                }
            }
            return *this;       
        }  
        
        inline mat_hgsw &get_PHI_matrix(int idx, int inv) { // copy contents of constant PHI matrix to the given matrix
            if(inv == 0){
                for (uint64_t i = 0; i < ROW; i++) {
                    for (uint64_t j = 0; j < COLUMN; j++){
                        this->mat[i][j] = prm::PHI_19[idx][i][j];
                    }
                }    
            }else{
                for (uint64_t i = 0; i < ROW; i++) {
                    for (uint64_t j = 0; j < COLUMN; j++){
                        this->mat[i][j] = prm::PHI_19_pr[idx][i][j];
                    }
                }
            }
            return *this;
        }

        template <uint32_t C2>
        inline mat_hgsw &mul_NTT_qpr(const mat_hgsw<T, prm, ROW, C2> &mat1, const mat_hgsw<T, prm, C2, COLUMN> &mat2) { // special function to perform the operation at Appx.F Alg.4 Line 5
            mat_hgsw<T, prm, ROW, COLUMN> res_first, res_second;
            mat_hgsw<T, prm, ROW, C2> mat1_first, mat1_second;
            mat_hgsw<T, prm, C2, COLUMN> mat2_first, mat2_second;
            for (uint32_t k = 0; k < ROW; k++) { 
                for (uint32_t i = 0; i < C2; i++) { 
                    for (uint32_t j = 0; j < prm::d; j++) {     
                        mat1_first[k][i][j] = mat1[k][i][j] % prm::p;  
                        mat1_second[k][i][j] = mat1[k][i][j] % prm::qprop; 
                    }
                    NTT(mat1_first[k][i], prm::qbar, prm::rp_conv_fact_qbar, prm::root_qbar);
                    NTT(mat1_second[k][i], prm::qprop, prm::rp_conv_fact_qprop, prm::root_qprop);
                }
            }        
            for (uint32_t k = 0; k < C2; k++) { 
                for (uint32_t i = 0; i < COLUMN; i++) { 
                    for (uint32_t j = 0; j < prm::d; j++) {     
                        mat2_first[k][i][j] = mat2[k][i][j] % prm::p;  
                        mat2_second[k][i][j] = mat2[k][i][j] % prm::qprop;     
                    }
                    NTT(mat2_first[k][i], prm::qbar, prm::rp_conv_fact_qbar, prm::root_qbar);
                    NTT(mat2_second[k][i], prm::qprop, prm::rp_conv_fact_qprop, prm::root_qprop);
                }
            }        
            res_first = mat1_first * mat2_first;
            res_second.mul_mod_qprop(mat1_second, mat2_second); 
            for (uint32_t k = 0; k < ROW; k++) {
                for (uint32_t i = 0; i < COLUMN; i++) {
                    INTT(res_first[k][i], prm::qbar, prm::inv_n_qbar, prm::root_inv_qbar);            
                    for (uint32_t j = 0; j < prm::d; j++) {
                        uint64_t sign = (res_first[k][i][j] >> 58) * 0x7ffffffffffffff;
                        res_first[k][i][j] = ((res_first[k][i][j] - (sign & prm::qbar)) + prm::q) % prm::p; 
                    } 
                    INTT(res_second[k][i], prm::qprop, prm::inv_n_qprop, prm::root_inv_qprop);
                }
                for (uint32_t i = 0; i < COLUMN; i++) { 
                    for (uint32_t j = 0; j < prm::d; j++) {
                        uint64_t tmp = ((res_first[k][i][j] - res_second[k][i][j]) + prm::q) % prm::p;
                        uint64_t h = (tmp * prm::qprop_inv) % prm::p;
                        this->mat[k][i][j] = res_second[k][i][j] + h * prm::qprop;
                    }
                }       
            }
            return *this;       
        }
          
    };

    template <typename Zq, typename zp, class prm>
    void sample_matrix_uniform(pair<mat_hgsw<Zq, prm, prm::mq, prm::n>, mat_hgsw<Zq, prm, prm::mq, prm::n>> *At, aes256ctr_ctx ctx, uint16_t t) { // fill the given matrix by performing uniform sampling, work in NTT form
        for (uint32_t i = 0; i < prm::m[t]; i++) {
            for (uint32_t j = 0; j < prm::mq; j++) {
                for (uint32_t k = 0; k < prm::n; k++) {
                    uni_sampler_zzp_1d_out64(At[i].first[j][k], prm::d, &ctx);
                    NTT(At[i].first[j][k], prm::qbar, prm::rp_conv_fact_qbar, prm::root_qbar);
                    uni_sampler_zzqbar_1d(At[i].second[j][k], prm::d, &ctx);  
                }
            }
        }
    }

    template <typename Zq, typename zp, class prm>
    void sample_matrix_gaussian(pair<mat_hgsw<Zq, prm, prm::mq, prm::ell_pr>, mat_hgsw<Zq, prm, prm::mq, prm::ell_pr>> *Et, DiscreteGaussianGenerator dgg, uint16_t t) { // fill the given matrix by performing discrete gaussian sampling, work in NTT form
        mat_hgsw<Zq, prm, prm::mq, prm::ell_pr> E[prm::m[t]];       
        for(uint64_t l = 0; l < prm::m[t]; l++) {
            for (uint32_t i = 0; i < prm::mq; i++) { 
                for (uint32_t j = 0; j < prm::ell_pr; j++) { 
                    for (uint32_t k = 0; k < prm::d; k++) {
                        E[l][i][j][k] = (dgg.GenerateIntegerKarney(0, prm::sig_s) + prm::q) % prm::q;      
                        Et[l].first[i][j][k] = E[l][i][j][k] % prm::p;    
                        E[l][i][j][k] = E[l][i][j][k] % prm::qbar;
                    }
                    NTT(Et[l].first[i][j], prm::qbar, prm::rp_conv_fact_qbar, prm::root_qbar);
                    NTT(E[l][i][j], prm::qbar, prm::rp_conv_fact_qbar, prm::root_qbar);
                    for (uint32_t k = 0; k < prm::d; k++) {     
                        Et[l].second[i][j][k] = E[l][i][j][k]; 
                    }            
                }
            }
        }
    }

    template <typename Zq, typename zp, class prm>
    void get_St(pair<mat_hgsw<Zq, prm, prm::n, prm::ell_pr>, mat_hgsw<Zq, prm, prm::n, prm::ell_pr>> &St, mat_hgsw<Zq, prm, prm::n, prm::ell_pr> &S) { // special function to create St in NTT form, using S generated in Setup
        mat_hgsw<Zq, prm, prm::n, prm::ell_pr> S_tmp;
        for (uint32_t k = 0; k < prm::n; k++) { 
            for (uint32_t i = 0; i < prm::ell_pr; i++) { 
                for (uint32_t j = 0; j < prm::d; j++) {     
                    S_tmp[k][i][j] = S[k][i][j];    
                }
            }
        }
        for (uint32_t k = 0; k < prm::n; k++) { 
            for (uint32_t i = 0; i < prm::ell_pr; i++) { 
                for (uint32_t j = 0; j < prm::d; j++) {     
                    St.first[k][i][j] = S_tmp[k][i][j] % prm::p;
                    S_tmp[k][i][j] = S_tmp[k][i][j] % prm::qbar;     
                }
                NTT(St.first[k][i], prm::qbar, prm::rp_conv_fact_qbar, prm::root_qbar);                
                NTT(S_tmp[k][i], prm::qbar, prm::rp_conv_fact_qbar, prm::root_qbar);
                for (uint32_t j = 0; j < prm::d; j++) {     
                    St.second[k][i][j] = S_tmp[k][i][j]; 
                }            
            }
        }
    }


    template <typename Zq, typename zp, class prm>
    auto NTT_mat(mat_hgsw<Zq, prm, 1, prm::mq> a_hatT_i_rand) { // special function to perform the operation at Appx.F Alg.3 Line 6 
        pair<mat_hgsw<Zq, prm, 1, prm::mq>, mat_hgsw<Zq, prm, 1, prm::mq>> atT_i_rand;    
        for (uint32_t i = 0; i < prm::mq; i++) { 
            for (uint32_t j = 0; j < prm::d; j++) {     
                atT_i_rand.first[0][i][j] = a_hatT_i_rand[0][i][j] % prm::p; 
                atT_i_rand.second[0][i][j] = a_hatT_i_rand[0][i][j] % prm::qbar;   
            }
            NTT(atT_i_rand.first[0][i], prm::qbar, prm::rp_conv_fact_qbar, prm::root_qbar);                        
            NTT(atT_i_rand.second[0][i], prm::qbar, prm::rp_conv_fact_qbar, prm::root_qbar);
        }
        return atT_i_rand;
    }

    template <typename Zq, class prm>
    auto INTT_mat(pair<mat_hgsw<Zq, prm, 1, prm::n + prm::ell_pr>, mat_hgsw<Zq, prm, 1, prm::n + prm::ell_pr>> ctstar) { 
        mat_hgsw<Zq, prm, 1, prm::n + prm::ell_pr> cstar; 
        for (uint32_t i = 0; i < prm::n + prm::ell_pr; i++) { 
            INTT(ctstar.first[0][i], prm::qbar, prm::inv_n_qbar, prm::root_inv_qbar);            
            for (uint32_t j = 0; j < prm::d; j++) {
                uint64_t sign = (ctstar.first[0][i][j] >> 58) * 0x7ffffffffffffff;
                ctstar.first[0][i][j] = ((ctstar.first[0][i][j] - (sign & prm::qbar)) + prm::q) % prm::p; 
            }            
            INTT(ctstar.second[0][i], prm::qbar, prm::inv_n_qbar, prm::root_inv_qbar);
            for (uint32_t j = 0; j < prm::d; j++) {
                uint64_t tmp = ((ctstar.first[0][i][j] - ctstar.second[0][i][j]) + prm::q) % prm::p;
                uint64_t h = (tmp * prm::qb_inv) % prm::p;
                cstar[0][i][j] = ctstar.second[0][i][j] + h * prm::qbar;
            }
        }
        return cstar;
    }

    void SampleD(double sigma, size_t k, double* c, double* d, size_t base, DiscreteGaussianGenerator dgg, int64_t* z) { // auxiliary function used in the gadget sampler
        z[k-1] = dgg.GenerateIntegerKarney(-c[k-1]/d[k-1],sigma/d[k-1]);
        for(size_t i = 0;i < k;++i){
            c[i] = c[i] - (static_cast<double>(z[k - 1])) * d[i];
        }
      for (size_t i = 0; i < k - 1; i++) {
        z[i] = dgg.GenerateIntegerKarney(-c[i], sigma);
      }    
    }

    void Perturb_normal_(double sigma, size_t k, double* l, double* h, double* p) { // auxiliary function used in the gadget sampler
        double z[k];
        std::normal_distribution<> d(0, sigma);

        PRNG &g = PseudoRandomNumberGenerator::GetPRNG();
       
        for(size_t i = 0;i < k;i++){
                z[i] = d(g);
        }
        for (size_t i = 0; i < k - 1; i++){
            p[i] = l[i]*z[i]+h[i+1]*z[i+1];        
        }
        p[k-1] = h[k-1]*z[k-1];
    }
    
    /* Implementation of the gadget sampler by Yang Yu: https://eprint.iacr.org/2021/1664.pdf */
    void GaussianABG_normal_(double stddev, size_t k, int64_t base, int64_t u, unsigned long long q, int64_t* t, DiscreteGaussianGenerator dgg) {
        double sigma = stddev / ((double)base + 1);
        int64_t q_digit[k];
        unsigned long long qq = q;
        for(size_t i = 0;i < k;++i){
            q_digit[i] = qq % base;
            qq = (unsigned long long)(qq / base);
        }
        int64_t u_digit[k];
        int64_t uu = u;
        for(size_t i = 0;i < k;++i){ 
            u_digit[i] = uu % base;
            uu = (int64_t)(uu / base);
        }
        double l[k],h[k];
        l[0] = sqrt(base * (1 + 1 / k) + 1);
        for(size_t i = 1;i < k;++i){
            l[i] = sqrt(base * (1 + 1 / static_cast<double>(k - i)));
        }
        
        h[0] = 0;
        for (size_t i = 1; i < k; i++)
            h[i] = sqrt(base * (1 - 1 / static_cast<double>(k - (i - 1))));
        double d[k];
        d[0] = ((int64_t)q_digit[0]) / static_cast<double>(base);
        for (size_t i = 1; i < k; i++){
            d[i] = (d[i-1] + (int64_t)q_digit[i]) / static_cast<double>(base);    
        }

        double p[k];
     
        Perturb_normal_(sigma,k,l,h,p);

        double c[k];
        c[0] = ((u_digit[0]) - p[0]) / static_cast<double>(base);
        for (size_t t = 1; t < k; t++) {
          c[t] = (c[t-1] + (u_digit[t]) - p[t]) /
                    static_cast<double>(base);
        }

        int64_t z[k];
        SampleD(sigma,k,c,d,base,dgg,z);

        t[0] = base * z[0] + (int64_t)(q_digit[0]) * z[k - 1] +
                     (int64_t)(u_digit[0]);
        t[k-1] = (int64_t)(q_digit[k-1]) * z[k-1] - z[k-2]+(int64_t)(u_digit[k-1]);
        for (size_t i = 1; i < k - 1; i++) {
          t[i] = base * z[i] - z[i - 1] +
                       (int64_t)(q_digit[i]) * z[k - 1] + (int64_t)(u_digit[i]);
        }  
    }  

    template <typename Zq, typename zp, class prm>
    auto g_inv_rand(mat_hgsw<zp, prm, prm::ell, 1> a) { // The g^inv_rand algorithm stated in Definition 3
        DiscreteGaussianGenerator dggKarney;
        mat_hgsw<Zq, prm, 1, prm::mq> z_mat;
        for(uint32_t i = 0; i < prm::d; i++){
            int64_t t[prm::mq];
            GaussianABG_normal_(prm::sig_r, prm::mq, prm::beta, a[0][0][i], prm::q, t, dggKarney);
            for(uint32_t j = 0; j < prm::mq; j++){ 
                z_mat[0][j][i] = (t[j] + prm::q) % prm::q;                     
            }        
        }
        return z_mat;
    }
    
    template <typename Zq, class prm>
    auto ModSwitch(mat_hgsw<Zq, prm, 1, prm::n + prm::ell_pr> &cstar) { // perform modulo switching from q to q', Appx.F Alg.3 Line 15
        mat_hgsw<Zq, prm, 1, prm::n + prm::ell_pr> cstarpr;
        for (uint32_t i = 0; i < prm::n + prm::ell_pr; i++) {
            for (uint32_t j = 0; j < prm::d; j++) { 
                double dblval = cstar[0][i][j];
                cstarpr[0][i][j] = (uint64_t)round(prm::qproq * dblval);
            }
        }
        return cstarpr;
    }    
    
    template <typename zp_t, class prm>
    auto Eval_phi(mat_hgsw<zp_t, prm, prm::f, prm::f> &PHI, zp_t *aX) { // auxiliary function to correctly perform CRT Encode
        mat_hgsw<zp_t[prm::d], prm, 1, 1> res;
        for(uint32_t i = 0; i < prm::f; i++){
            uint32_t mtmp[prm::f];
            for(uint32_t j = 0; j < prm::f; j++){
                 mtmp[j] = (aX[i] * PHI[i][j]) % prm::p;
            }
            for(uint32_t j = 0; j < prm::f; j++){
                res[0][0][j] = (res[0][0][j] + mtmp[j]) % prm::p;
            }    
        }
        return res;
    }

    template <typename zp_t, class prm>
    auto CRT_poly(mat_hgsw<zp_t[prm::d], prm, 1, 1> *residues, mat_hgsw<zp_t[prm::d], prm, 1, 4> &f_arr) { // auxiliary function to correctly perform CRT Encode
        mat_hgsw<zp_t[prm::d], prm, 1, 1> res;
        zp_t del[prm::d];
        for(uint32_t j = 0; j < prm::f; j++){
            del[j] = ((residues[0][0][0][j] - residues[1][0][0][j]) + prm::p) % prm::p;
        }
        for(uint32_t j = prm::f; j < prm::d; j++){
            del[j] = 0;
        }
        zp_t c1[prm::d], c2[prm::d];
        rp_mulp(c1, del, f_arr[0][2]); 
        rp_mulp(c2, c1, f_arr[0][0]);
        for(uint32_t j = 0; j < prm::d; j++){
            res[0][0][j] = ((residues[0][0][0][j] - c2[j]) + prm::p) % prm::p;
        }
        return res;
    }     
    
    template <typename zp_t, uint32_t size, class prm>
    auto CRT_Encode(mat_hgsw<zp_t[prm::d], prm, 1, prm::w> &qi, mat_hgsw<zp_t, prm, prm::f, prm::f> *PHI, mat_hgsw<zp_t[prm::d], prm, 1, 4> &f_arr) { // CRT encoding used both in Encrypt and Add
        mat_hgsw<zp_t[prm::d], prm, 1, 1> r_hat[prm::ell];
        uint32_t v = size / prm::ell_p;
        for(uint32_t j = 0; j < v; j++){
            mat_hgsw<zp_t[prm::d], prm, 1, 1> residues[prm::ell_p];
            for(uint32_t i = 0; i < prm::ell_p; i++){
                residues[i] = Eval_phi<zp_t, prm>(PHI[i], qi[0][j * prm::ell_p + i]);
            }
            r_hat[j] = CRT_poly<zp_t, prm>(residues, f_arr);
            
        }
        mat_hgsw<zp_t[prm::d], prm, prm::ell, 1> res;
        for(uint32_t j = 0; j < prm::ell; j++){
            for(uint32_t k = 0; k < prm::d; k++){
                res[j][0][k] = r_hat[j][0][0][k];       
            }
        }
        return res;
    }  

    template <typename zp_t, uint32_t v, class prm>
    auto CRT_Decode(mat_hgsw<zp_t[prm::d], prm, prm::ell, 1> &mu_bar, mat_hgsw<zp_t, prm, prm::f, prm::f> *PHI_pr, mat_hgsw<zp_t[prm::d], prm, 1, 4> &f_arr) { // CRT decoding used both in Decrypt
        mat_hgsw<zp_t[prm::d], prm, 1, 1> r_eval[v * prm::ell_p];
        mat_hgsw<zp_t[prm::d], prm, v, prm::ell_p> res;

        for(uint32_t j = 0; j < v; j++) {
            for(uint32_t i = 0; i < prm::ell_p; i++) {
                uint32_t deg_r = 0;    
                for(int32_t l = prm::d-1; l >= 0; l--) {
                    if(mu_bar[j][0][l] != 0) {
                        deg_r = l;
                        break;
                    }        
                }
                zp_t r[prm::d] = {0};
                for(uint32_t l = 0; l < prm::d; l++) {
                    r[l] = mu_bar[j][0][l];
                }
                uint32_t s_val = 0, s_loc = 0;
                while(deg_r >= prm::f) {
                    for(int32_t l = prm::d-1; l >= 0; l--) {
                        if(r[l] != 0) {
                            s_val = r[l];   
                            s_loc = deg_r - prm::f;  
                            break;
                        }        
                    }
                    zp_t sb[prm::d] = {0};
                    for(uint32_t l = 0; l < prm::f+1; l++) {
                        if(f_arr[0][i][l] != 0) {
                            sb[s_loc + l] = (s_val * f_arr[0][i][l]) % prm::p;
                        }
                    }
                    for(uint32_t l = 0; l < prm::d; l++) {
                        r[l] = ((r[l] - sb[l]) + prm::p) % prm::p;
                    }
                    for(int32_t l = prm::d-1; l >= 0; l--) {
                        if(r[l] != 0) {
                            deg_r = l;
                            break;
                        }
                    }
                }        
                r_eval[j * prm::ell_p + i] = Eval_phi<zp_t, prm>(PHI_pr[i], r);    
            }
        }   
        for(uint32_t i = 0; i < v; i++) {
            for(uint32_t j = 0; j < prm::ell_p; j++) {
                for(uint32_t k = 0; k < prm::f; k++) {
                    res[i][j][k] = r_eval[i * prm::ell_p + j][0][0][k];
                }
            }
        }
        return res;        
    }
    
    /* validate the result of Decrypt(q_res) by comparing it with a_i * q_i, for i = 1, 2, ..., w */
    template <typename zp, class prm>
    void validate(mat_hgsw<zp, prm, 1, 1> *a, mat_hgsw<zp, prm, 1, prm::w> *q_, mat_hgsw<zp, prm, prm::v, prm::ell_p> &q_res, mat_hgsw<zp, prm, 1, 4> &f_arr, uint16_t t) {
        mat_hgsw<zp, prm, prm::v, prm::ell_p> aq_;
        for(uint32_t i = 0; i < prm::m[t]; i++) {
            for(uint32_t j = 0; j < prm::v; j++) {
                for(uint32_t k = 0; k < prm::ell_p; k++) {   
                    zp r = {0};
                    rp_mulp(r, a[i][0][0], q_[i][0][j * prm::ell_p + k]);       
                    uint32_t deg_r = 0;
                    for(int32_t l = prm::d-1; l >= 0; l--) {
                        if(r[l] != 0) {
                            deg_r = l;
                            break;
                        }        
                    }
                    uint32_t s_val = 0, s_loc = 0;
                    while(deg_r >= prm::f) {
                        for(int32_t l = prm::d-1; l >= 0; l--) {
                            if(r[l] != 0) {
                                s_val = r[l];   
                                s_loc = deg_r - prm::f;  
                                break;
                            }        
                        }
                        zp sb = {0};
                        for(uint32_t l = 0; l < prm::f+1; l++) {
                            if(f_arr[0][0][l] != 0) {
                                sb[s_loc + l] = (s_val * f_arr[0][0][l]) % prm::p;
                            }
                        }
                        for(uint32_t l = 0; l < prm::d; l++) {
                            r[l] = ((r[l] - sb[l]) + prm::p) % prm::p;
                        }
                        for(int32_t l = prm::d-1; l >= 0; l--) {
                            if(r[l] != 0) {
                                deg_r = l;
                                break;
                            }
                        }
                    }
                    for(uint32_t l = 0; l < prm::f; l++) {
                        aq_[j][k][l] = (aq_[j][k][l] + r[l]) % prm::p;
                    }
                }
            }            
        }
        uint32_t fl = 1;
        for (uint32_t i = 0; i < prm::v; i++) {
            for (uint32_t j = 0; j < prm::ell_p; j++) {
                for (uint32_t l = 0; l < prm::f; l++) {
                    if(q_res[i][j][l] != aq_[i][j][l]) {
                        fl = 0;
                    }
                }   
            }
        }
        if(fl == 1) {
            cout << "\nPASS\n";
        }else {
            cout << "\nFAIL\n";
            q_res.print_matrix();
            printf("\n\n");
            aq_.print_matrix();
        }
    }
}
#endif

