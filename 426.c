#include <stdio.h>
#include <pbc.h>
#include "prime.h"
#include "pbc_functions_inshort.h"
void main(){
	int number_of_iot = 400;
	number_of_iot -=1;

	double start,end;

	init_pairing();

	t g,h;
	t msk;
	IG1(g);
	IG1(h);
	Zr(msk);
	//element_printf("%B\n",g);
	//element_printf("%B\n",h);

	/*	
		double start, end;
		Zr(h);
		start = pbc_get_time();
		for(int i=0;i<100;i++) 	element_from_hash(g,h,100);
		end = pbc_get_time();
		printf(" Z->G time : %fms\n",end-start);

		start = pbc_get_time();
		for(int i=0;i<100;i++) element_from_hash(h,g,100);
		end = pbc_get_time();
		printf(" G->Z time : %fms\n",end-start);

		t gt;
		IGT(gt);
		start = pbc_get_time();
		for(int i=0;i<100;i++) element_from_hash(h,gt,100);
		end = pbc_get_time();
		printf(" GT->Z time : %fms\n",end-start);
		*/
/******************/
	t x,Y;
	Zr(x);
	IG1(Y);
	mpz_t N,N_square,p_,q_,f1,f2,x0,lambda,mu;
	mpz_init(N);
	mpz_init(p_);
	mpz_init(q_);
	mpz_init(f1);
	mpz_init(f2);
	mpz_init(x0);
	mpz_init(lambda);
	mpz_init(mu); 
t EID;
Zr(EID);
t sk0;
IG1(sk0);
	t HID[10000];
	t ski[10000];
	t J_HIDj[10000];
	for(int i=0;i<=number_of_iot;i++)
	{
		Zr(HID[i]);
		IG1(J_HIDj[i]);
		IG1(ski[i]);
	}
	t Ts,sigma_cs;
	Zr(Ts);
	IG1(sigma_cs);

	t H_Y_EID_Ts;
	IG1(H_Y_EID_Ts);
	t temp1,temp2,temp3;
	IG1(temp1);
	IG1(temp2);
	IG1(temp3);
	t left,right;	
	IGT(left);
	IGT(right);
t mi[10000];
	t alphai[10000];
	t betai[10000];
	t Ci[10000];
	t Ri[10000];
	IG1(temp1);
	IG1(temp2);
	IG1(temp3);
	t H_Ts;
	IG1(H_Ts);
	for(int i=0;i<=number_of_iot;i++)
	{	
		Zr(mi[i]);

		//4.4.1
		Zr(alphai[i]);
		Zr(betai[i]);
		IG1(Ci[i]);
		IG1(Ri[i]);
	}
t Kij[4000][10];
	t psi_Kij[4000][10];
	

	t si[10000];
	for(int i=0;i<=number_of_iot;i++){
		IGT(Kij[i][0]);
		IGT(Kij[i][1]);
		Zr(psi_Kij[i][0]);
		Zr(psi_Kij[i][1]);
		Zr(si[i]);
	}
	t Vi[1000];
	mpz_t Di[1000];
	mpz_t Ri_[1000];
	mpz_t mi_mpz[1000];
	mpz_t betai_mpz[1000];
	for(int i=0;i<=number_of_iot;i++){
		mpz_init(Di[i]);
		mpz_init(mi_mpz[i]);
		mpz_init(betai_mpz[i]);
		mpz_init(Ri_[i]);
		_mpz_random(Ri_[i],256);	//256 bits exponents
		IG1(Vi[i]);
	}
t s0,K0j[2],psi_K0j[2];
	Zr(s0);
	//	element_set0(s0);
	for(int j=0;j<=1;j++){
		IGT(K0j[j]);
		Zr(psi_K0j[j]);
	}
	t CVi[1000];
	for(int i=0;i<=number_of_iot;i++){
		IG1(CVi[i]);
	}
t h0_i[1000];
	mpz_t h0_i_[1000];
	mpz_t Ri__[1000];
	for(int i=0;i<=number_of_iot;i++){
		mpz_init(h0_i_[i]);
		Zr(h0_i[i]);
		mpz_init(Ri__[i]);
	}
	mpz_t p_1,q_1,p_1_q_1;
	mpz_init(p_1);
	mpz_init(q_1);
	mpz_init(p_1_q_1);
	//4.6.2
	
	for(int i=0;i<=number_of_iot;i++){
		mpz_init(h0_i_[i]);
		Zr(h0_i[i]);
		mpz_init(Ri__[i]);
	}

	start=pbc_get_time();

	
	/************************************************************/
	double start_part,end_part,start_test,end_test;
// 4.2 Registration
	// 4.2.1
	start_part=pbc_get_time();
	element_pow_zn(Y,h,x); 
	mpz_t *point_p_ = &p_;
	mpz_t *point_q_ = &q_;
	prime_generate(point_p_,point_q_);
	mpz_mul(N,p_,q_);
	mpz_pow_ui(N_square,N,2);
	mpz_set(mu,p_);
	//mpz_set(x0,mu);
	_mpz_random(x0,256);
	mpz_powm_ui(f1,mu,2,N_square);
	mpz_powm(f2,f1,x0,N_square);
	
	//gmp_printf("p_ is : %Zd\n",p_);
	//gmp_printf("q_ is : %Zd\n",q_);
	   

	//4.2.2
	element_from_hash(sk0,EID,100);	
	element_pow_zn(sk0,sk0,msk);

	//4.2.3
	for(int i=0;i<=number_of_iot;i++)
	{
		element_from_hash(J_HIDj[i],HID[i],100);
		element_pow_zn(ski[i],J_HIDj[i],msk);
	}
	end_part=pbc_get_time();
	printf("registration:%.3f\n",(end_part-start_part)*1000);
// 4.3 Data Collection Request
	start_part=pbc_get_time();
	element_from_hash(temp1,Y,100);
	element_from_hash(temp2,EID,100);
	element_from_hash(temp3,Ts,100);
	element_add(H_Y_EID_Ts,temp1,temp2);
	element_add(H_Y_EID_Ts,H_Y_EID_Ts,temp3);
	element_clear(temp1);
	element_clear(temp3);
	element_clear(temp2);

	//element_printf("H_Y_EID_Ts:%B\n",H_Y_EID_Ts);//
	//element_printf("x: %B\n",x);//
	element_pow_zn(sigma_cs,H_Y_EID_Ts,x);
	//element_printf("sigma:%B\n",sigma_cs);//
	//4.3.2
	element_pairing(left,sigma_cs,h);
	//element_printf("%B\n",h);//
	//element_printf("%B\n",left);//
	element_pairing(right,H_Y_EID_Ts,Y);
	//element_printf("%B\n",right);//
//	if(element_cmp(left,right)==0) printf("4.3.2 holds!\n");

	end_part=pbc_get_time();
	printf("data cloolection request:%.3f\n",(end_part-start_part)*1000);
	//4.4 Hybrid IoT Devices Report Generation
	start_part=pbc_get_time();
	element_from_hash(H_Ts,Ts,100);
	/*******
	for(int i=0;i<=number_of_iot;i++)
	{	
		element_printf("m%d=%B\n",i,mi[i]);
		element_printf("alphai%d=%B\n",i,alphai[i]);
		element_printf("betai%d=%B\n",i,betai[i]);
		element_printf("Ci%d=%B\n",i,Ci[i]);
	}
	*******/
	IG1(temp1);
	IG1(temp2);
	IG1(temp3);
	for(int i=0;i<=number_of_iot;i++)
	{
		element_pow_zn(temp1,g,mi[i]);
		element_pow_zn(temp2,H_Ts,alphai[i]);
		element_pow_zn(temp3,Y,betai[i]);
		element_mul(Ci[i],temp1,temp2);
		element_mul(Ci[i],Ci[i],temp3);
		element_pow_zn(Ri[i],h,betai[i]);
		
	}
	element_clear(temp1);
	element_clear(temp2);
	element_clear(temp3);


	// 4.4.4
	
	
	Zr(temp1);
	Zr(temp2);
	for(int i=0;i<=number_of_iot;i++){
		element_set0(temp1);
		element_set0(temp2);
		//IGT(Kij[i][0]);
		//IGT(Kij[i][1]);
		//Zr(psi_Kij[i][0]);
		//Zr(psi_Kij[i][1]);
		//for(int j=0;j<=number_of_iot;j++){
		//	if(j<i-1) continue;
		//	if(j>i+1) break;

		//printf("[%d,%d]\n",i,j);//
		//if(j!=i){
		//printf("[%d,%d]\n",i,j);//
		if(i==0) {
			element_pairing(Kij[i][0],J_HIDj[number_of_iot],ski[i]);
			element_pairing(Kij[i][1],J_HIDj[i+1],ski[i]);	
		}
		if(i==number_of_iot) {
			element_pairing(Kij[i][0],J_HIDj[i-1],ski[i]);
			element_pairing(Kij[i][1],J_HIDj[0],ski[i]);
		}
		if(i>0&&i<number_of_iot){
			element_pairing(Kij[i][0],J_HIDj[i-1],ski[i]);
			element_pairing(Kij[i][1],J_HIDj[i+1],ski[i]);
		}
		//printf("[%d,%d]\n",i,j);//
		element_from_hash(psi_Kij[i][0],Kij[i][0],100);
		element_from_hash(psi_Kij[i][1],Kij[i][1],100);
		//Zr(si[i]);
		element_sub(si[i],psi_Kij[i][1],psi_Kij[i][0]);
		//element_printf("s%d = %B\n",i,si[i]);//
		//element_printf("%B\n",psi_Kij[i][j]);//
		/*
		   if(j<i){
		   element_add(temp1,temp1,psi_Kij[i][j]);
		   }
		   if(j>i){
		   element_add(temp2,temp2,psi_Kij[i][j]);
		   }
		   */
		//	}

	
	}
	/*
	t s0test,s0test_;
	Zr(s0test);
	Zr(s0test_);
	element_set0(s0test);
	for(int i=0;i<number_of_iot;i++) element_add(s0test,si[i],s0test);
			element_printf("s0set = %B\n",s0test);//element_sub(s0test,s0test,s0test
			element_neg(s0test_,s0test);
			element_add(s0test,s0test,s0test_);
			element_printf("s0set-s0set = %B\n",s0test);//
			*/
	//4.4.5

	mpz_t temp4;
	mpz_init(temp4);
	Zr(temp1);
	for(int i=0;i<=number_of_iot;i++){
		//mpz_init(Di[i]);
		//mpz_init(mi_mpz[i]);
		//mpz_init(betai_mpz[i]);
		//mpz_init(Ri_[i]);
		//_mpz_random(Ri_[i],256);	//256 bits exponents
		//IG1(Vi[i]);
		element_to_mpz(mi_mpz[i],mi[i]);
		//element_to_mpz(betai_mpz[i],betai[i]);
		_mpz_random(betai_mpz[i],256);
		mpz_mul(temp4,N,mi_mpz[i]);
		mpz_add_ui(temp4,temp4,1);
		mpz_powm(Di[i],f2,betai_mpz[i],N_square);	
		mpz_mul(Di[i],Di[i],temp4);

		element_from_hash(temp1,Ci[i],100); //没有计算量 
		element_to_mpz(Ri_[i],temp1);

		mpz_mul(Ri_[i],Ri_[i],betai_mpz[i]);
		mpz_powm(Ri_[i],f1,Ri_[i],N_square);

		element_sub(temp1,si[i],alphai[i]);
		element_pow_zn(Vi[i],H_Ts,temp1);
	}
	/**************************************************************
	2020-05-01 17:37:06 

	**************************************************************/
	end_part=pbc_get_time();
	printf("Hybrid IoT Devices Report Generation:%.3f\n",(end_part-start_part)*1000);
//4.5
	start_part=pbc_get_time();
//4.5.1

//4.5.2
//compute s0
	
	//	for(int j=0;j<=number_of_iot;j++) {
	//	}
	//	for(int j=0;j<=number_of_iot;j++){
	element_pairing(K0j[0],sk0,J_HIDj[number_of_iot-1]);
	element_pairing(K0j[1],sk0,J_HIDj[1]);
	element_from_hash(psi_K0j[0],K0j[0],100);
	element_from_hash(psi_K0j[1],K0j[1],100);
	element_sub(s0,psi_K0j[0],psi_K0j[1]);
/*
	   t a;
	   Zr(a);

	   element_printf("s0=%B\n",a);
	   element_neg(a,a);	
	   element_printf("s0=%B\n",a);

*/
	t C0;
	IG1(C0);
	element_pow_zn(C0,H_Ts,s0);
	//4.5.3
	t L1;
	IG1(L1);
	element_clear(temp1);
	IG1(temp1);
	element_set1(temp1);
	for(int i=0;i<=number_of_iot;i++){
//		IG1(CVi[i]);
		element_mul(CVi[i],Ci[i],Vi[i]);
		element_mul(temp1,temp1,CVi[i]);
	}
	element_mul(L1,C0,temp1);

	t L2;
	IG1(L2);
	element_set1(L2);
	for(int j=0;j<=number_of_iot;j++) element_mul(L2,L2,Ri[j]);

	mpz_t D;
	mpz_init(D);
	mpz_set_ui(D,1);
	for(int j=0;j<=number_of_iot;j++) mpz_mul(D,D,Di[j]);
	end_part=pbc_get_time();
	printf("4.5:%.3f\n",(end_part-start_part)*1000);

//4.6
	start_part=pbc_get_time();
	mpz_sub_ui(p_1,p_,1);	
	mpz_sub_ui(q_1,q_,1);	
	mpz_mul(p_1_q_1,p_1,q_1);
	for(int i=0;i<=number_of_iot;i++){
		//mpz_init(h0_i_[i]);
		//Zr(h0_i[i]);
		//mpz_init(Ri__[i]);
		element_from_hash(h0_i[i],CVi[i],100);
		element_to_mpz(h0_i_[i],h0_i[i]);
		mpz_invert(h0_i_[i],h0_i_[i],p_1_q_1);
		mpz_powm(Ri__[i],Ri_[i],h0_i_[i],N_square);
	}	
	//4.6.3
	mpz_t	L3;
	mpz_init(L3);
	mpz_set_ui(L3,1);
	for(int i=0;i<=number_of_iot;i++){
		mpz_mul(L3,L3,Ri__[i]);
	}

	mpz_t M;
	mpz_init(M);
	mpz_powm(M,L3,x0,N_square);
	mpz_cdiv_q(M,D,M);
	mpz_mod(M,M,N_square);
	mpz_sub_ui(M,M,1);
	mpz_cdiv_q(M,M,N); 
	mpz_mod(M,M,N_square);
	// up to now, the decryption process is completed
	element_pow_mpz(temp1,g,M);
	element_pow_zn(temp2,L2,x);
	element_mul(temp1,temp1,temp2);
	//element_printf("left = %B\n",temp1);
	//element_printf("right = %B\n",L1);

	/*****************
	mpz_t n2time,n2time_exp;
	mpz_init(n2time);
	_mpz_random(n2time,512);
	mpz_init(n2time_exp);
	_mpz_random(n2time_exp,256);
	t ptime;
	Zr(ptime);
	double start,end;
	start=pbc_get_time();
	for(int i=1;i<=60000;i++) mpz_powm(n2time,n2time,n2time_exp,N_square);
	end=pbc_get_time();
	printf("Z_N^2 exponent time: %f\n",(end-start)/60);
	start=pbc_get_time();
	for(int i=1;i<=60000;i++) element_pow_mpz(ptime,ptime,q_);
	end=pbc_get_time();
	printf("Z_p exponent time: %f\n",(end-start)/60);

	*****************/
	end_part=pbc_get_time();
	printf("4.6:%.3f\n",(end_part-start_part)*1000);
	end=pbc_get_time();
	printf("number:%d\ntotal time=%.3f\n",number_of_iot+1,(end-start)*1000);
	}
