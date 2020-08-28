#include <stdio.h>
#include <pbc.h>
#include "prime.h"
#include "pbc_functions_inshort.h"
void main(){
	int number_of_iot = 1000;
	number_of_iot -=1;


	init_pairing();

	t g,h;
	t msk;
	IG1(g);
	IG1(h);
	Zr(msk);


double test1,test2;
/********************/
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
	_mpz_random(x0,512);
	_mpz_random(p_,512);
	mpz_nextprime(p_,p_);
	mpz_nextprime(q_,p_);

	t EID;
	Zr(EID);
	t sk0;
	IG1(sk0);

	t HID[10000];
	t ski[10000];
	t J_HIDj[10000];
	for(int j=0;j<=number_of_iot;j++){
		Zr(HID[j]);
		IG1(J_HIDj[j]);
		IG1(ski[j]);
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
	mpz_t temp4;
	mpz_init(temp4);
	
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
		mpz_init(Di[i]);
		mpz_init(mi_mpz[i]);
		mpz_init(betai_mpz[i]);
		mpz_init(Ri_[i]);
		_mpz_random(Ri_[i],160);	//256 bits exponents
		_mpz_random(betai_mpz[i],512);


		Zr(alphai[i]);
		Zr(betai[i]);
		IG1(Ci[i]);
		IG1(Ri[i]);
	}

	t s0,K0j[2],psi_K0j[2];
	Zr(s0);
	for(int j=0;j<=1;j++){
		IGT(K0j[j]);
		Zr(psi_K0j[j]);
	}


	t C0;
	IG1(C0);
	
	t L1;
	IG1(L1);

	t L2;
	IG1(L2);

	mpz_t D;
	mpz_init(D);

	mpz_t h0,R_;
	mpz_init_set_ui(R_,1);
	mpz_init_set_ui(h0,1);

	mpz_t temp5;
	mpz_init_set_ui(temp5,1);
	mpz_t H0_Ci[10000];
	for(int i=0;i<=number_of_iot;i++){
		mpz_init(H0_Ci[i]);
		_mpz_random(H0_Ci[i],160);
	}

	mpz_t h0_invert;
	mpz_t M;
	mpz_init(h0_invert);
	mpz_init(M);

	double iot1,iot2,iot3,node1,node2,node3,server1,server2,server3,total1,total2,total3,agg1,agg2,agg3,de1,de2,de3;
/*********************/
// 4.2 Registration
	// 4.2.1
	total1 = pbc_get_time();
	server1 = pbc_get_time();
	element_pow_zn(Y,h,x); 
	
	mpz_mul(N,p_,q_);
	mpz_pow_ui(N_square,N,2);
	mpz_set(mu,p_);
	mpz_powm_ui(f1,mu,2,N_square);
	mpz_powm(f2,f1,x0,N_square);
	server2 = pbc_get_time();
	server3 = server2-server1;	
	   

	//4.2.2
	node1 = pbc_get_time();
	element_from_hash(sk0,EID,100);	
	element_pow_zn(sk0,sk0,msk);
	node2 = pbc_get_time();
	node3 = node2-node1;

	//4.2.3
	iot1 = pbc_get_time();
	for(int j=0;j<=number_of_iot;j++){
		element_from_hash(J_HIDj[j],HID[j],100);
		element_pow_zn(ski[j],J_HIDj[j],msk);
	}
	iot2 = pbc_get_time();
	iot3 = iot2 - iot1;
// 4.3 Data Collection Request
	
	server1= pbc_get_time();
	element_from_hash(H_Y_EID_Ts,Y,100);

	element_clear(temp1);
	element_clear(temp3);
	element_clear(temp2);

	element_pow_zn(sigma_cs,H_Y_EID_Ts,x);
	server2= pbc_get_time();
	server3 = server2 - server1;
	//4.3.2
	node1 = pbc_get_time();
	element_pairing(left,sigma_cs,h);
	element_from_hash(H_Y_EID_Ts,Y,100);
	element_pairing(right,H_Y_EID_Ts,Y);
	node2 = pbc_get_time();
	node3 = node2-node1;
	//if(element_cmp(left,right)==0) printf("4.3.2 holds!\n");

//6.1 Hybrid IoT Devices Report Generation
	//6.1.3
iot1 = pbc_get_time();	

	for(int i=0;i<=number_of_iot;i++){
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
		element_from_hash(psi_Kij[i][0],Kij[i][0],100);
		element_from_hash(psi_Kij[i][1],Kij[i][1],100);
		element_sub(si[i],psi_Kij[i][1],psi_Kij[i][0]);
	}


	//6.1.4
	element_from_hash(H_Ts,Ts,100);

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
	for(int i=0;i<=number_of_iot;i++){
		element_to_mpz(mi_mpz[i],mi[i]);
		mpz_mul(temp4,N,mi_mpz[i]);
		mpz_add_ui(temp4,temp4,1);
		mpz_powm(Di[i],f2,betai_mpz[i],N_square);	
		mpz_mul(Di[i],Di[i],temp4);


		mpz_mul(Ri_[i],Ri_[i],betai_mpz[i]);
		mpz_powm(Ri_[i],f1,Ri_[i],N_square);

	}
	element_clear(temp1);
	element_clear(temp2);
	element_clear(temp3);
	iot2 = pbc_get_time();
	iot3 = iot3 + iot2 - iot1;

//6.2
//compute s0
	agg1 = pbc_get_time();
	node1 = pbc_get_time();
	element_pairing(K0j[0],sk0,J_HIDj[number_of_iot-1]);
	element_pairing(K0j[1],sk0,J_HIDj[1]);
	element_from_hash(psi_K0j[0],K0j[0],100);
	element_from_hash(psi_K0j[1],K0j[1],100);
	element_sub(s0,psi_K0j[0],psi_K0j[1]);


	element_pow_zn(C0,H_Ts,s0);
	//6.2.2
	element_set(L1,C0);
	for(int i=0;i<=number_of_iot;i++){
		element_mul(L1,Ci[i],L1);
	}
	element_set1(L2);
	for(int j=0;j<=number_of_iot;j++) element_mul(L2,L2,Ri[j]);

	mpz_set_ui(D,1);
	for(int j=0;j<=number_of_iot;j++) mpz_mul(D,D,Di[j]);

	Zr(temp1);	
       /////////////
       test1 = pbc_get_time();	
	for(int i=0;i<=number_of_iot;i++){
		mpz_mul(h0,h0,H0_Ci[i]);
		mpz_mod(h0,h0,N_square);
	}
    
	test2 = pbc_get_time();
	printf("h0:%f\n",(test2-test1)*1000);
       test1 = pbc_get_time();	
	for(int i=0;i<=number_of_iot;i++){
		mpz_cdiv_q(temp5,h0,H0_Ci[i]);
	} 
	test2 = pbc_get_time();
	printf("h0:%f\n",(test2-test1)*1000);
	//test1 = pbc_get_time();/////////////////
	mpz_t temp_mpz;
	mpz_init(temp_mpz);
	for(int i=0;i<=number_of_iot;i++){
		mpz_powm(temp_mpz,Ri_[i],temp5,N_square);
	}
	//test2 = pbc_get_time();////////////////
	//printf("test:%.3f\n\n\n",(test2-test1)*1000); /////////////////
	for(int i=0;i<=number_of_iot;i++){
		mpz_mul(R_,R_,temp5);
		mpz_mod(R_,R_,N_square);
	}	


	agg2 = pbc_get_time();
	agg3 = agg3 +agg2 - agg1;
	node2 = pbc_get_time();
	node3 = node3 +node2 - node1;
	
	//6.3Decryption
	element_clear(temp1);
	IG1(temp1);
	IG1(temp2);
	de1 = pbc_get_time();
	server1 = pbc_get_time();
	mpz_invert(h0_invert,h0,N_square);	
	//mpz_mod(h0_invert,h0_invert,N_square);	
	mpz_mul(temp4,x0,h0_invert);
	//mpz_mod(R_,R_,N_square);	
	//mpz_mod(temp4,temp4,N_square);	
	mpz_powm(temp4,R_,temp4,N_square);
	test1 = pbc_get_time();
	mpz_cdiv_q(temp4,D,temp4);
	mpz_mod(temp4,temp4,N_square);
	test2 = pbc_get_time();
	printf("test2-test1:%f\n",(test2-test1)*1000);
	mpz_cdiv_q(M,temp4,N);
	mpz_mod(M,M,N_square);


	element_pow_mpz(temp1,g,M);
	element_pow_zn(temp2,L2,x);
	element_mul(temp1,temp1,temp2);
	server2 = pbc_get_time();
	server3 = server3+server2-server1;

	de2 = pbc_get_time();
	de3 = de2-de1;

	total2 = pbc_get_time();
	total3 = total2 - total1;
	printf("number:%d\n",number_of_iot+1);
	printf("total:%.3f\n",total3*1000);
	printf("iot:%.3f\n",iot3*1000/number_of_iot);
	printf("node:%.3f\n",node3*1000);
	printf("server:%.3f\n",server3*1000);
	printf("agg:%.3f\n",agg3*1000);
	printf("de:%.3f\n",de3*1000);
	
/*	
	
double start,end;
_mpz_random(N,256);

start = pbc_get_time();
for(int i=0;i<6000;i++){
	mpz_cdiv_q(p_,temp4,N);
	//mpz_cdiv_q(p_,q_,N,N_square);
}
end = pbc_get_time();
printf("divn2:%.3f\n",(end-start)/6);

_mpz_random(N,1024);

start = pbc_get_time();
for(int i=0;i<6000;i++){
	mpz_powm(p_,N,N,N_square);
}
end = pbc_get_time();
printf("1024En2:%.3f\n",(end-start)/6);
_mpz_random(N,512);

start = pbc_get_time();
for(int i=0;i<6000;i++){
	mpz_powm(p_,N,N,N_square);
}
end = pbc_get_time();
printf("512En2:%.3f\n",(end-start)/6);
//Mn2
start = pbc_get_time();
for(int i=0;i<60000;i++){
mpz_mul(p_,N,N);
mpz_mod(p_,p_,N_square);
}
end = pbc_get_time();
printf("Mn2:%.3f\n",(end-start)/60);
//Eg1
element_clear(temp1);
element_clear(temp2);
IG1(temp1);
Zr(temp2);
start = pbc_get_time();
for(int i=0;i<6000;i++){
	element_pow_zn(temp1,temp1,temp2);	
}
end = pbc_get_time();
printf("Eg1:%.3f\n",(end-start)/6);
//Mg1
start = pbc_get_time();
for(int i=0;i<6000;i++){
element_mul(temp1,temp1,temp1);	
}
end = pbc_get_time();
printf("Mg1:%.3f\n",(end-start)/6);
//Egt
start = pbc_get_time();
for(int i=0;i<6000;i++){
	element_pow_zn(left,left,temp2);	
}
end = pbc_get_time();
printf("Egt:%.3f\n",(end-start)/6);
//P
test1=0.0;
element_random(temp1);
for(int i=0;i<6000;i++){
start = pbc_get_time();
	element_pairing(left,temp1,temp1);	
end = pbc_get_time();
test1+=end-start;
}
printf("p:%.3f\n",(test1)/6);
//Hm2p
start = pbc_get_time();
for(int i=0;i<6000;i++){
	element_from_hash(temp1,temp2,100);	
}
end = pbc_get_time();
printf("H2p:%.3f\n",(end-start)/6);
*/

}
