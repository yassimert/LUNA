#ifndef __HGSW_ALGOS__
#define __HGSW_ALGOS__

#include "lattice/dgsampling.h"
#include "hgsw_funcs.hpp"

using namespace hgsw_params;
using namespace lbcrypto;
using namespace hgsw_funcs;

using namespace std::chrono;
namespace hgsw_algos {

   template <typename Zq, typename zp, class Params>
   void HGSW_Setup( mat_hgsw<Zq, Params, Params::n, Params::ell_pr> &S,
                    pair<mat_hgsw<Zq, Params, Params::n, Params::ell_pr>, mat_hgsw<Zq, Params, Params::n, Params::ell_pr>> &St,
                    mat_hgsw<zp, Params, Params::tau, Params::ell> &T,
                    pair<mat_hgsw<Zq, Params, Params::mq, Params::n>, mat_hgsw<Zq, Params, Params::mq, Params::n>> *At,
                    pair<mat_hgsw<Zq, Params, Params::mq, Params::ell_pr>, mat_hgsw<Zq, Params, Params::mq, Params::ell_pr>> *Et, uint16_t t) {
               
	    uint8_t seed[32];
	    aes256ctr_ctx ctx;
	    randombytes(seed, Params::d);
	    aes256ctr_init(&ctx, seed, 0);
        DiscreteGaussianGenerator dgg(Params::s);
        
        sample_matrix_uniform<Zq, zp, Params>(At, ctx, t);
        sample_matrix_gaussian<Zq, zp, Params>(Et, dgg, t);
       
        T.sample_matrix_uniform_p();
        S.sample_matrix_gaussian(dgg);
        get_St<Zq, zp, Params>(St, S);        
    }

    template <typename Zq, typename zp, class Params>
    auto HGSW_Encrypt( pair<mat_hgsw<Zq, Params, Params::n, Params::ell_pr>, mat_hgsw<Zq, Params, Params::n, Params::ell_pr>> &St,
                       mat_hgsw<zp, Params, Params::tau, Params::ell> &T,
                       pair<mat_hgsw<Zq, Params, Params::mq, Params::n>, mat_hgsw<Zq, Params, Params::mq, Params::n>> &Ait,
                       pair<mat_hgsw<Zq, Params, Params::mq, Params::ell_pr>, mat_hgsw<Zq, Params, Params::mq, Params::ell_pr>> &Eit,    
                       mat_hgsw<zp, Params, 1, Params::w> &qi,
                       mat_hgsw<uint16_t, Params, Params::f, Params::f> *PHI,
                       mat_hgsw<zp, Params, 1, 4> &f_arr,                       
                       mat_hgsw<Zq, Params, 1, Params::mq> &gT_p ) {

        mat_hgsw<zp, Params, Params::ell, 1> mui = CRT_Encode<uint16_t, Params::w, Params>(qi, PHI, f_arr);    
        mat_hgsw<zp, Params, 1, 1> Tmui; 
        Tmui.mul_modp(T, mui);
        
        mat_hgsw<zp, Params, 1, Params::ell_pr> mu_bar_i;
        
        for(uint32_t i = 0; i < Params::ell; i++) {
            for(uint32_t j = 0; j < Params::d; j++) { 
                mu_bar_i[0][i][j] = mui[0][i][j];
            }
        }  
        for(uint32_t j = 0; j < Params::d; j++) { 
            mu_bar_i[0][Params::ell][j] = Tmui[0][0][j];
        }
        
        pair<mat_hgsw<Zq, Params, Params::mq, Params::ell_pr>, mat_hgsw<Zq, Params, Params::mq, Params::ell_pr>> AitSt_Eit;

        AitSt_Eit.first = Ait.first * St.first;
        AitSt_Eit.second = Ait.second * St.second;

        AitSt_Eit.first += Eit.first; 
        AitSt_Eit.second += Eit.second; 
        
        mat_hgsw<Zq, Params, 1, Params::ell_pr> mu_bar_ntt;
        for(uint32_t i = 0; i < Params::ell_pr; i++) {
            for (uint32_t k = 0; k < Params::d; k++) {
                mu_bar_ntt[0][i][k] = mu_bar_i[0][i][k];     
            }
            NTT(mu_bar_ntt[0][i], Params::qbar, prm::rp_conv_fact_qbar, prm::root_qbar);
        }
        
        mat_hgsw<Zq, Params, Params::mq, Params::ell_pr> Hit_p;        
        Hit_p.compute_Hit_p(mu_bar_ntt, gT_p);
        AitSt_Eit.first += Hit_p;
         
        pair<mat_hgsw<Zq, Params, Params::mq, Params::n + Params::ell_pr>, mat_hgsw<Zq, Params, Params::mq, Params::n + Params::ell_pr>> Cit; 
        for(uint32_t i = 0; i < Params::mq; i++) {
            for(uint32_t j = 0; j < Params::n; j++) {
                for(uint32_t l = 0; l < Params::d; l++) {
                    Cit.first[i][j][l] = Ait.first[i][j][l];
                    Cit.second[i][j][l] = Ait.second[i][j][l];
                }
            }
            uint32_t k = 0;
            for(uint32_t j = Params::n; j < Params::n + Params::ell_pr; j++) {
                for(uint32_t l = 0; l < Params::d; l++) {
                    Cit.first[i][j][l] = AitSt_Eit.first[i][k][l]; 
                    Cit.second[i][j][l] = AitSt_Eit.second[i][k][l]; 
                 }
                 k += 1;        
            }
        }
        return Cit;
    }

    template <typename Zq, typename zp, class Params>
    auto HGSW_Add( pair<mat_hgsw<Zq, Params, Params::mq, Params::n + Params::ell_pr>, mat_hgsw<Zq, Params, Params::mq, Params::n + Params::ell_pr>> *Ct,
                   mat_hgsw<zp, Params, 1, 1> *a,
                   mat_hgsw<uint16_t, Params, Params::f, Params::f> *PHI,
                   mat_hgsw<zp, Params, 1, 4> &f_arr, uint16_t t) {
        
        DiscreteGaussianGenerator dgg;
        pair<mat_hgsw<Zq, Params, 1, Params::n + Params::ell_pr>, mat_hgsw<Zq, Params, 1, Params::n + Params::ell_pr>> ctstar;
        for(uint32_t j = 0; j < Params::m[t] / Params::nu; j++) {       
            for(uint32_t i = 0; i < Params::nu; i++) { 
                mat_hgsw<zp, Params, 1, Params::w> a_arr; 
                for(uint32_t k = 0; k < Params::ell_p; k++) {
                    for(uint32_t l = 0; l < Params::f; l++) {
                        a_arr[0][k][l] = a[j * Params::nu + i][0][0][l];
                    }
                }
                mat_hgsw<zp, Params, Params::ell, 1> a_hat_i = CRT_Encode<uint16_t, Params::ell_p, Params>(a_arr, PHI, f_arr);    

                mat_hgsw<Zq, Params, 1, Params::mq> a_hatT_i_rand;
                a_hatT_i_rand = g_inv_rand<Zq, zp, Params>(a_hat_i);

                pair<mat_hgsw<Zq, Params, 1, Params::mq>, mat_hgsw<Zq, Params, 1, Params::mq>> atT_i_rand;
                atT_i_rand = NTT_mat<Zq, zp, Params>(a_hatT_i_rand);

                pair<mat_hgsw<Zq, Params, 1, Params::n + Params::ell_pr>, mat_hgsw<Zq, Params, 1, Params::n + Params::ell_pr>> tmp;            
                tmp.first = atT_i_rand.first * Ct[j * Params::nu + i].first;
                tmp.second = atT_i_rand.second * Ct[j * Params::nu + i].second; 
                
                ctstar.first += tmp.first;
                ctstar.second += tmp.second;
            }
        }
        mat_hgsw<Zq, Params, 1, Params::n + Params::ell_pr> cstar; 
        cstar = INTT_mat<Zq, Params>(ctstar);

        for(uint32_t j = 0; j < Params::m[t] / Params::nu; j++) {
            mat_hgsw<Zq, Params, 1, Params::n + Params::ell_pr> zero_n_yT;
            for(uint32_t i = Params::n; i < Params::n + Params::ell_pr; i++) {  
                for (uint32_t k = 0; k < Params::d; k++) {
                    zero_n_yT[0][i][k] = (dgg.GenerateIntegerKarney(0, Params::sig_r) + Params::q) % Params::q;
                }
            }
            cstar.add_mod(zero_n_yT, Params::q);
        }
        
        return ModSwitch<Zq, Params>(cstar);      
    }

    template <typename Zq, typename zp, class Params>
    auto HGSW_Decrypt( mat_hgsw<Zq, Params, Params::n, Params::ell_pr> &S,
                       mat_hgsw<zp, Params, Params::tau, Params::ell> &T,
                       mat_hgsw<Zq, Params, 1, Params::n + Params::ell_pr> &cstarpr,
                       mat_hgsw<uint16_t, Params, Params::f, Params::f> *PHI_pr,
                       mat_hgsw<zp, Params, 1, 4> &f_arr ) {

        mat_hgsw<Zq, Params, Params::ell_pr, Params::n + Params::ell_pr> S_bar;        
        for(uint32_t i = 0; i < Params::ell_pr; i++) {  
            for(uint32_t j = 0; j < Params::n; j++) {  
                for(uint32_t k = 0; k < Params::d; k++) {  
                    S_bar[i][j][k] = (-S[j][i][k] + Params::q) % Params::q;
                }
            }
        }
        mat_hgsw<Zq, Params, Params::ell_pr, Params::ell_pr> I_ell_pr;  
        for(uint32_t i = 0; i < Params::ell_pr; i++) {  
            for(uint32_t j = 0; j < Params::ell_pr; j++) {  
                if(i == j){
                    I_ell_pr[i][j][0] = 1;
                }else{
                    I_ell_pr[i][j][0] = 0; 
                }
                for(uint32_t k = 1; k < Params::d; k++) {  
                    I_ell_pr[i][j][k] = 0;
                } 
            }
        }
        for(uint32_t i = 0; i < Params::ell_pr; i++) {  
            uint32_t k = 0;
            for(uint32_t j = Params::n; j < Params::n + Params::ell_pr; j++) {  
                for(uint32_t l = 0; l < Params::d; l++) {     
                    S_bar[i][j][l] = I_ell_pr[i][k][l]; 
                }
                k += 1;
            }
        }        
        mat_hgsw<Zq, Params, Params::n + Params::ell_pr, Params::ell_pr> S_barT;        
        for(uint32_t i = 0; i < Params::n + Params::ell_pr; i++) {  
            for(uint32_t j = 0; j < Params::ell_pr; j++) {  
                for(uint32_t l = 0; l < Params::d; l++) {   
                    S_barT[i][j][l] = S_bar[j][i][l];
                }
            }
        } 
        
        mat_hgsw<Zq, Params, 1, Params::ell_pr> H_bar, H_barzp;  
        
        for (uint32_t i = 0; i < Params::n + Params::ell_pr; i++) {
            for (uint32_t j = 0; j < Params::ell_pr; j++) {
                for (uint32_t l = 0; l < Params::d; l++) {
                    uint64_t sign = (S_barT[i][j][l] >> 62) * 0x7fffffffffffffff;
                    S_barT[i][j][l] = ((S_barT[i][j][l] - (sign & Params::q)) + Params::qpr) % Params::qpr; 
                } 
            }
        }
        
        H_bar.mul_NTT_qpr(cstarpr, S_barT);
        H_barzp.scale_convert(H_bar);
        
        mat_hgsw<zp, Params, Params::ell, 1> mu_bar_1;  
        mat_hgsw<zp, Params, 1, Params::tau> mu_bar_2;  
        for(uint32_t i = 0; i < Params::ell; i++) {    
            for(uint32_t j = 0; j < Params::d; j++) {
                mu_bar_1[i][0][j] = H_barzp[0][i][j];   
            }       
        }
        for(uint32_t j = 0; j < Params::d; j++) {
            mu_bar_2[0][0][j] = H_barzp[0][Params::ell][j];   
        } 
         
        mat_hgsw<zp, Params, 1, Params::tau> Tmu_bar_1;  
        Tmu_bar_1.mul_modp(T, mu_bar_1);
       
        for(uint32_t l = 0; l < Params::d; l++) {
            if(Tmu_bar_1[0][0][l] != mu_bar_2[0][0][l]) {
                cout << "\nFAIL\n";
                Tmu_bar_1.print_matrix();         
                mu_bar_2.print_matrix(); 
                exit(1);
            }
        }

        return CRT_Decode<uint16_t, Params::v, Params>(mu_bar_1, PHI_pr, f_arr); 
    }  

}
#endif
