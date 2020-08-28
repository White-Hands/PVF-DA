#include <stdio.h>
#include "pbc_functions_inshort.h"
#include "prime.h"
void main(){
	double start,end;
	double start_part,end_part;
	double start_test,end_test;
	double iot1,iot2,iot3,node1,node2,node3,server1,server2,server3,agg1,agg2,agg3,de1,de2,de3,total1,total2,total3;
	int number_of_iot;
	number_of_iot=1000;
	mpz_t p,q,p_sub_1,q_sub_1,n,n_square,lambda,mu;
	mpz_init(p);
	mpz_init(p_sub_1);
	mpz_init(q);
	mpz_init(q_sub_1);
	mpz_init(n);
	mpz_init(n_square);
	mpz_init(lambda);
	mpz_init(mu);

	_mpz_random(p,512);
	mpz_nextprime(p,p);
	mpz_nextprime(q,p);
	mpz_mul(n,p,q);
	mpz_sub_ui(p_sub_1,p,1);	
	mpz_sub_ui(q_sub_1,q,1);	
	mpz_mul(lambda,p_sub_1,q_sub_1);
	mpz_pow_ui(n_square,n,2);
	mpz_t g;
	mpz_init(g);
	_mpz_random(g,512);
	mpz_powm(mu,g,lambda,n_square);	
	mpz_sub_ui(mu,mu,1);
	mpz_cdiv_q(mu,mu,n);
	mpz_invert(mu,mu,n_square);

//	gmp_printf("mu is %Zd\n",mu);

	init_pairing();
	t alpha_hat,x_hat,g1;
	Zr(alpha_hat);	
	Zr(x_hat);	
	IG1(g1);

	t Q;
	IG1(Q);
	t g1_g1_a_hat;
	IGT(g1_g1_a_hat);
	element_pairing(g1_g1_a_hat,g1,g1);
	element_pow_zn(g1_g1_a_hat,g1_g1_a_hat,alpha_hat);

	t Y_hat;
	IG1(Y_hat);
	element_pow_zn(Y_hat,g1,x_hat);
	t omega;
	Zr(omega);
//参数初始化
	t up,down;
	IGT(up);
	IGT(down);
	t x;
	Zr(x);
	t Omega_Res;
	IG1(Omega_Res);
	t C2_hat,C3_hat;
	IG1(C3_hat);
	IG1(C2_hat);
	t C1_hat;
	IGT(C1_hat);

	t MR,beta_hat;
	IGT(MR);
	Zr(beta_hat);

t y,z;
	Zr(y);
	Zr(z);
	t g2,g3;
	IG1(g2);
	IG1(g3);
	t si[number_of_iot+1];
	t ui[number_of_iot+1];
	t Hchi[number_of_iot+1];
	t BLS[number_of_iot+1];
	for(int i=0;i<number_of_iot;i++){
		Zr(si[i]);
		Zr(ui[i]);
	IG1(Hchi[i]);
	IG1(BLS[i]);
	}
	t Xi[number_of_iot+1];
	t Yi[number_of_iot+1];
	t ki[number_of_iot+1];
	t ri[number_of_iot+1];
	t alphai[number_of_iot+1];
	t betai[number_of_iot+1];
	for(int i=0;i<number_of_iot;i++){
		Zr(Xi[i]);
		IG1(Yi[i]);
		Zr(ki[i]);
		Zr(ri[i]);

		IG1(alphai[i]);
		Zr(betai[i]);
	}
	t alpha_yanzhengi[number_of_iot+1];
	t temp1;
	IG1(temp1);
	for(int i=0;i<number_of_iot;i++){
		IG1(alpha_yanzhengi[i]);
	}
	t ti[number_of_iot+1];
	t aki1[number_of_iot+1];
	t aki2[number_of_iot+1];
	t aki3[number_of_iot+1];
	for(int i=0;i<number_of_iot;i++){
		Zr(ti[i]);
		IG1(aki1[i]);
		IG1(aki2[i]);
		IG1(aki3[i]);
	}
	t temp2;
	IGT(temp2);
	t temp3;
	IGT(temp3);
	mpz_t ci[number_of_iot+1];
	mpz_t mi[number_of_iot+1];
	mpz_t vi[number_of_iot+1];
	mpz_t temp_mpz_1;
	mpz_init(temp_mpz_1);
	//start_test=pbc_get_time();
	for(int i=0;i<number_of_iot;i++){
	mpz_init(ci[i]);
	mpz_init(mi[i]);
	mpz_init(vi[i]);
	_mpz_random(vi[i],512);
	_mpz_random(mi[i],512);
	}
	t si_[number_of_iot];
	t ui_[number_of_iot];
	t z_invert;
	Zr(z_invert);
	t temp_zn;
	Zr(temp_zn);
	for(int i=0;i<number_of_iot;i++){
		Zr(si_[i]);
		Zr(ui_[i]);
	}
	t Xj[number_of_iot];
	t Omega_Agg[number_of_iot];
	for(int i=0;i<1;i++){
		Zr(Xj[i]);
		IG1(Omega_Agg[i]);
	}



	start = pbc_get_time();
//B Registration
//User Registration
	start_part = pbc_get_time();
	iot1=pbc_get_time();
	for(int i=0;i<number_of_iot;i++){
		element_pow_zn(Yi[i],g1,Xi[i]);
		element_pow_zn(alphai[i],g1,ri[i]);
		element_mul(betai[i],Xi[i],Xi[i]);
		element_sub(betai[i],ri[i],betai[i]);
	}
	iot2=pbc_get_time();
	iot3 = pbc_get_time();
// Authentication
//printf("Authentication\n");
	server1 = pbc_get_time();

	for(int i=0;i<number_of_iot;i++){
		element_pow_zn(alpha_yanzhengi[i],g1,betai[i]);
		element_pow_zn(temp1,Yi[i],betai[i]);
		element_mul(alpha_yanzhengi[i],alpha_yanzhengi[i],temp1);
	}

	for(int i=0;i<number_of_iot;i++){
		element_pow_zn(aki1[i],g1,alpha_hat);
		element_pow_zn(temp1,Yi[i],ti[i]);
		element_mul(aki1[i],aki1[i],temp1);

		element_pow_zn(aki2[i],Q,ti[i]);
		element_pow_zn(aki3[i],g1,ti[i]);
	}
	server2 = pbc_get_time();
	server3 = server2 - server1;
//Offline Signature Generation
	iot1 = pbc_get_time();
	for(int i=0;i<number_of_iot;i++){
	element_pow_zn(g2,g1,y);
	element_pow_zn(g3,g1,z);
	element_pow_zn(temp1,g1,ri[i]);

	element_pow_zn(Hchi[i],g2,si[i]);
	element_mul(Hchi[i],Hchi[i],temp1);
	element_pow_zn(temp1,g3,ui[i]);
	element_mul(Hchi[i],Hchi[i],temp1);
	element_from_hash(Hchi[i],Hchi[i],100);

	element_pow_zn(BLS[i],Hchi[i],Xi[i]);
	}

	iot2 = pbc_get_time();
	iot3 = iot2-iot1+iot3;
	end_part=pbc_get_time();
//	printf("registration:%.3f\n",(end_part-start_part)*1000);
//Report Generation
start_part =pbc_get_time();

	node1=pbc_get_time();

	for(int i=0;i<number_of_iot;i++){
		element_from_hash(Hchi[i],Hchi[i],100);
		element_pairing(temp2,Yi[i],Hchi[i]);
		element_mul(temp3,temp3,temp2);
	}

	for(int i=0;i<number_of_iot;i++){
		element_mul(temp1,temp1,BLS[i]);
	}
	element_pairing(temp2,g1,temp1);

		node2=pbc_get_time();
		node3 = node2-node1;

		iot1 = pbc_get_time();
	for(int i=0;i<number_of_iot;i++){
		mpz_powm(ci[i],g,mi[i],n_square);
		mpz_powm(temp_mpz_1,vi[i],n,n_square);
		mpz_mul(ci[i],ci[i],temp_mpz_1);
		mpz_mod(ci[i],ci[i],n_square);
	}
	

	element_invert(z_invert,z);
	for(int i=0;i<number_of_iot;i++){
	//	Zr(si_[i]);
	//	Zr(ui_[i]);
		element_sub(ri[i],ri[i],ui[i]);
		element_sub(temp_zn,si[i],si_[i]);
		element_mul(temp_zn,temp_zn,y);
		element_mul(temp_zn,ui[i],z);
		element_add(temp_zn,ri[i],temp_zn);
		element_add(temp_zn,temp_zn,temp_zn);
		element_pow_zn(ui_[i],temp_zn,z_invert);
	}
	end_part=pbc_get_time();
	//printf("report generation:%.3f\n",(end_part-start_part)*1000);
		iot2 = pbc_get_time();
		iot3 = iot2-iot1+iot3;

//D Report Aggregation
start_part = pbc_get_time();
	node1 = pbc_get_time();
	for(int i=0;i<number_of_iot;i++){
		element_from_hash(temp1,ci[i],100);
	}
//Report Aggregation
//printf("//Report Aggregation\n");
	
	mpz_t c;
	mpz_init_set_ui(c,1);
	for(int i=0;i<number_of_iot;i++){
		mpz_mul(c,c,ci[i]);
		mpz_mod(c,c,n_square);
	}
//Aggregation Signature Generation
//printf("//Aggregation Signature Generation\n");

	for(int i=0;i<1;i++){
		//Zr(Xj[i]);
		//IG1(Omega_Agg[i]);
		element_from_hash(Omega_Agg[i],c,100);
		element_pow_zn(Omega_Agg[i],Omega_Agg[i],Xj[i]);
	}
	end_part=pbc_get_time();
	agg3 = end_part-start_part;
	//printf("report aggregation:%.3f\n",(end_part-start_part)*1000);
	node2 = pbc_get_time();
	node3 = node3+node2-node1;
//Report Reading
start_part = pbc_get_time();
de1 = pbc_get_time();
////Aggregation Signature Verification
	for(int i=0;i<1;i++){
		element_pairing(temp2,g1,Omega_Agg[i]);

		element_from_hash(Omega_Agg[i],c,100);
		element_pow_zn(Omega_Agg[i],Omega_Agg[i],Xj[i]);
		element_pairing(temp2,g1,Omega_Agg[i]);
	}
// Report Reading and Decryption
	for(int i=0;i<1;i++){
		mpz_powm(temp_mpz_1,c,lambda,n_square);
		mpz_sub_ui(temp_mpz_1,temp_mpz_1,1);
		mpz_cdiv_q(temp_mpz_1,temp_mpz_1,n);
		mpz_powm(temp_mpz_1,g,lambda,n_square);
		mpz_sub_ui(temp_mpz_1,temp_mpz_1,1);
		mpz_cdiv_q(temp_mpz_1,temp_mpz_1,n);
	}

	end_part=pbc_get_time();
	server3 = server3+end_part-start_part;
	de2 = pbc_get_time();
	de3 = de2 -de1;

//Response
	start_part=pbc_get_time();
	server1 = pbc_get_time();
//Step-1
	element_pow_zn(C1_hat,C1_hat,beta_hat);
	element_mul(C1_hat,C1_hat,MR);

	element_pow_zn(C2_hat,g1,beta_hat);
	element_div(Q,Yi[1],Q);
	element_pow_zn(C2_hat,Q,beta_hat);

	element_from_hash(Omega_Res,C1_hat,100);
	element_pow_zn(Omega_Res,Omega_Res,x);	
	server2 = pbc_get_time();
	server3 = server2-server1+server3;
//Step-2
	node1= pbc_get_time();
	element_pairing(temp2,g1,Omega_Res);
	element_from_hash(temp1,C2_hat,100);
	element_pairing(temp2,temp1,C2_hat);
	node2= pbc_get_time();
	node3 = node3+node2-node1;
//Step-3
//printf("//Step-3\n");
	iot1 = pbc_get_time();
	for(int i=0;i<number_of_iot;i++){
	element_pairing(up,C2_hat,aki1[i]);
	element_pairing(temp2,C2_hat,aki2[i]);
	element_pairing(down,C3_hat,aki3[i]);
	element_mul(down,temp2,temp2);
	element_div(temp2,up,down);

	element_div(temp2,C1_hat,temp2);

	}
	iot2= pbc_get_time();
	iot3 = iot3+iot2-iot1;
	end_part=pbc_get_time();
	printf("response:%.3f\n",(end_part-start_part)*1000);
	
	end = pbc_get_time();
	printf("number:%d\n",number_of_iot);
	printf("total:%.3f\n",(end-start)*1000);
	printf("iot:%.3f\n",iot3*1000/number_of_iot);
	printf("node:%.3f\n",node3*1000);
	printf("server:%.3f\n",server3*1000);
	printf("agg:%.3f\n",agg3*1000);
	printf("de:%.3f\n",de3*1000);
	
}



