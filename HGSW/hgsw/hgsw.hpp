#ifndef __HGSW__
#define __HGSW__

#include <chrono>
#include "hgsw_algos.hpp"

using namespace std::chrono;
using namespace hgsw_algos;
namespace hgsw {

    template <class Params>
    void HGSW() { // Half-GSW Encryption Scheme
        for(uint16_t t = 0; t < 3; t++) { // Ng = 2^10, 2^12, 2^14
            cout << "\nRunning with Ng = " << Params::m[t]/2 << "\n\n";

            mat_hgsw<uint16_t, Params, Params::f, Params::f> PHI[Params::f_num], PHI_pr[Params::f_num];
            for(uint16_t i = 0; i < Params::f_num; i++) { 
                PHI[i].get_PHI_matrix(i, 0); // get phi matrix used in CRT encode
                PHI_pr[i].get_PHI_matrix(i, 1); // get inverse phi matrix used in CRT decode
            }
            
            mat_hgsw<uint16_t[Params::d], Params, 1, 4> f_arr;
            f_arr.get_f_arr(); // get polynomials used in CRT encode
            
            mat_hgsw<uint64_t[Params::d], Params, 1, Params::mq> gT_p;        
            gT_p.gen_gT_p(); // get gT vector used in ginv_rand function

            mat_hgsw<uint64_t[Params::d], Params, Params::n, Params::ell_pr> S, S_tmp;

            pair<mat_hgsw<uint64_t[Params::d], Params, Params::n, Params::ell_pr>,
                 mat_hgsw<uint64_t[Params::d], Params, Params::n, Params::ell_pr>> St;
     
            mat_hgsw<uint16_t[Params::d], Params, Params::tau, Params::ell> T;
     
            pair<mat_hgsw<uint64_t[Params::d], Params, Params::mq, Params::n>,
                 mat_hgsw<uint64_t[Params::d], Params, Params::mq, Params::n>> At[Params::m[t]];
            
            pair<mat_hgsw<uint64_t[Params::d], Params, Params::mq, Params::ell_pr>,
                 mat_hgsw<uint64_t[Params::d], Params, Params::mq, Params::ell_pr>> Et[Params::m[t]];          
     
            auto start_setup = high_resolution_clock::now(); 
     
            HGSW_Setup<uint64_t[Params::d], uint16_t[Params::d], Params>(S, St, T, At, Et, t); // call Setup to create secret key sk = (S, T, A, E)

            auto end_setup = high_resolution_clock::now();
            auto time_setup = duration<double>(end_setup - start_setup);
            cout << "HGSW.Setup:" << time_setup.count() << " s" << "\n\n";

            mat_hgsw<uint16_t[Params::d], Params, 1, Params::w> q_[Params::m[t]];
            mat_hgsw<uint16_t[Params::d], Params, 1, 1> a[Params::m[t]];

            for(uint64_t i = 0; i < Params::m[t]; i++) {
                q_[i].gen_qi(); // generate plaintext vector qi
                a[i].gen_a(); // generate scaling factor a used in Add 
            }

            auto start_encrypt = high_resolution_clock::now();  
            
            pair<mat_hgsw<uint64_t[Params::d], Params, Params::mq, Params::n + Params::ell_pr>,
                 mat_hgsw<uint64_t[Params::d], Params, Params::mq, Params::n + Params::ell_pr>> Ct[Params::m[t]]; 

            for(uint64_t i = 0; i < Params::m[t]; i++) {
                Ct[i] = HGSW_Encrypt<uint64_t[Params::d], uint16_t[Params::d], Params>(St, T, At[i], Et[i], q_[i], PHI, f_arr, gT_p); // call Encrypt to compute ciphertexts Ct_i
            }

            auto end_encrypt = high_resolution_clock::now();          
            auto time_encrypt = duration<double>(end_encrypt - start_encrypt);
            cout << "HGSW.Encrypt:" << time_encrypt.count() << " s" << "\n\n";
      
            auto start_add = high_resolution_clock::now();  
         
            mat_hgsw<uint64_t[Params::d], Params, 1, Params::n + Params::ell_pr> cstarpr;
            cstarpr = HGSW_Add<uint64_t[Params::d], uint16_t[Params::d], Params>(Ct, a, PHI, f_arr, t);  // call Add to compute the re-randomized scaled ciphertext c*'     
           
            auto end_add = high_resolution_clock::now();          
            auto time_add = duration<double>(end_add - start_add);
            cout << "HGSW.Add:" << time_add.count() << " s" << "\n\n";
            
            auto start_decrypt = high_resolution_clock::now();  
           
            mat_hgsw<uint16_t[Params::d], Params, Params::v, Params::ell_p> q_res;       
            q_res = HGSW_Decrypt<uint64_t[Params::d], uint16_t[Params::d], Params>(S, T, cstarpr, PHI_pr, f_arr); // call Decrypt to recover the encrypted message  
            
            auto end_decrypt = high_resolution_clock::now();          
            auto time_decrypt = duration<double>(end_decrypt - start_decrypt);
            cout << "HGSW.Decrypt:" << time_decrypt.count() << " s" << "\n\n";

            validate<uint16_t[Params::d], Params>(a, q_, q_res, f_arr, t); // validate the result by comparing it with a_i * q_i, for i = 1, 2, ..., w
        }
    }

}
#endif

