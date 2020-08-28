#include <gmp.h>
#include <stdio.h>
#include <time.h>

void _mpz_random(mpz_t a,mp_bitcnt_t n)
{
    unsigned long int b;
    b = (unsigned long int)(time(NULL));
    //printf("time : %ld\n",b);
    //printf("Hello world!\n");
   // mpz_t a;
    //mpz_init(a);
    
    gmp_randstate_t grt;
    gmp_randinit_default(grt);
    //gmp_randseed_ui(grt,time(NULL));
    gmp_randseed_ui(grt,b);
    mpz_rrandomb(a,grt,n);
    }
    

void prime_generate(mpz_t *p,mpz_t *q){
    //unsigned long int b;
    //b = (unsigned long int)(time(NULL));
    //printf("time : %ld\n",b);
    //printf("Hello world!\n");
    mpz_t a;
    mpz_init(a);
    _mpz_random(a,512);
    
    //gmp_randstate_t grt;
    //gmp_randinit_default(grt);
    //gmp_randseed_ui(grt,time(NULL));
    //gmp_randseed_ui(grt,b);
    //mpz_rrandomb(a,grt,512);
    //gmp_printf ("%s is an mpz %Zd\n", "here", a);
        
    /*mpz_t p,q;
    mpz_init(p);
    mpz_init(q);
    */
    mpz_nextprime(*p,a);
    mpz_nextprime(*q,*p);
    //gmp_printf ("p is an mpz %Zd\n", p);    
    //gmp_printf ("q is %Zd\n",q);
    
    /*
    mpz_t n;
    mpz_init(n);
    mpz_mul(n,p,q);
    gmp_printf ("n is %Zd\n",n);
    mpz_t lambda,p_1,q_1;       
    mpz_init(lambda);
    mpz_init(q_1);
    mpz_init(p_1);
    mpz_sub_ui(p_1,p,1);    
    mpz_sub_ui(q_1,q,1);    
    mpz_mul(lambda,p_1,q_1); 
    gmp_printf ("lambda is %Zd\n",lambda);   */
}
