#include <stdio.h>
#include "pbc_functions_inshort.h"
#include "prime.h"
void main(){
	int number_of_iot;
	number_of_iot=1000;
	int number_of_type;
	number_of_type=1;
	double start,end,startpart,endpart,starttest,endtest;
	double iot1,iot2,iot3,total1,total2,total3,agg1,agg2,agg3,de1,de2,de3;

	init_pairing();
	mpz_t p1,q1,p_sub_1,q_sub_1,n,n_square,lambda,mu;
	mpz_init(p1);
	mpz_init(p_sub_1);
	mpz_init(q1);
	mpz_init(q_sub_1);
	mpz_init(n);
	mpz_init(n_square);
	mpz_init(lambda);
	mpz_init(mu);

	_mpz_random(p1,512);
	mpz_nextprime(p1,p1);
	mpz_nextprime(q1,p1);
	mpz_mul(n,p1,q1);
	mpz_sub_ui(p_sub_1,p1,1);	
	mpz_sub_ui(q_sub_1,q1,1);	
	mpz_mul(lambda,p_sub_1,q_sub_1);
	mpz_pow_ui(n_square,n,2);
	mpz_t g;
	mpz_init(g);
	_mpz_random(g,512);
	mpz_powm(mu,g,lambda,n_square);	
	mpz_sub_ui(mu,mu,1);
	mpz_cdiv_q(mu,mu,n);
	mpz_invert(mu,mu,n_square);

	/**************/
	t xg;
	Zr(xg);
	t Yg;
	IG1(Yg);
	t xi[number_of_iot];
	t Yi[number_of_iot];
	t ti1[number_of_iot];
	t ti2[number_of_iot];
	t aki1[number_of_iot];
	t aki2[number_of_iot];
	t aki3[number_of_iot];
	t aki4[number_of_iot];
	for(int i=0;i<number_of_iot;i++){
		Zr(xi[i]);
		Zr(ti1[i]);
		Zr(ti2[i]);
		IG1(Yi[i]);
		IG1(aki1[i]);
		IG1(aki2[i]);
		IG1(aki3[i]);
		IG1(aki4[i]);
	}
	mpz_t di[number_of_iot];
	mpz_t ri[number_of_iot];
	mpz_t Ci[number_of_iot];
	mpz_t temp_mpz;
	mpz_init(temp_mpz);
	for(int i=0;i<number_of_iot;i++){
		mpz_init(di[i]);
		mpz_init(ri[i]);
		_mpz_random(di[i],512);
		_mpz_random(ri[i],512);
		mpz_init_set_ui(Ci[i],1);
	}
	mpz_t Xj[number_of_type];	
	mpz_t Dj[number_of_type];	
	mpz_init(Xj[number_of_type-1]);
	for(int j=number_of_type-1;j>=2;j--){
		mpz_init(Xj[j-1]);
		mpz_init(Dj[j]);
	}
	t C1,C2,C3,C4;
	IGT(C1);
	IG1(C2);
	IG1(C3);
	IG1(C4);
	t m;
	IGT(m);
	t s;
	Zr(s);

	t up,down;
	IGT(up);
	IGT(down);

	t sigma;
	IG1(sigma);


	t temp_IG1;
	t temp_IGT;
	IG1(temp_IG1);
	IGT(temp_IGT);
	t P,Q1,Q2,Y;
	IG1(P);
	IG1(Q1);
	IG1(Q2);
	IG1(Y);

	t alpha,x;
	Zr(alpha);	
	Zr(x);	

	t ePPalpha;
	IGT(ePPalpha);

	double server1,server2,server3,user1,user2,user3,node1,node2,node3;
	/**************/
	//gmp_printf("mu is %Zd\n",mu);

// chao di zeng xu lie 
	server1=pbc_get_time();
	mpz_t ai[number_of_type];
	mpz_init(ai[0]);
	mpz_set_ui(ai[0],1);
	mpz_init(ai[1]);
	_mpz_random(ai[1],512);
	mpz_t gi[number_of_type];
	for(int i=0;i<number_of_type;i++){
		mpz_init(gi[i]);
		if(i>=2) {
			mpz_init(ai[i]);
			for(int j=0;j<i;j++){
				mpz_add(ai[i],ai[j],ai[i]);
			}
			mpz_nextprime(ai[i],ai[i]);
		}
		mpz_powm(gi[i],g,ai[i],n_square);
	}

	element_mul(Y,x,P);

	element_pairing(ePPalpha,P,P);
	element_pow_zn(ePPalpha,ePPalpha,alpha);
	server2 = pbc_get_time();
	server3=server2-server1;

// 开始！
// local gateway register
	node1 =pbc_get_time();
	element_pow_zn(Yg,P,xg);
	node2 =pbc_get_time();
	node3 = node2-node1;
// user register
	user1 = pbc_get_time();
	for(int i=0;i<number_of_iot;i++){
		//Zr(xi[i]);
		//Zr(ti1[i]);
		//Zr(ti2[i]);
		//IG1(Yi[i]);
		//IG1(aki1[i]);
		//IG1(aki2[i]);
		//IG1(aki3[i]);
		IG1(aki4[i]);
		element_pow_zn(Yi[i],P,xi[i]);
		element_mul(aki1[i],alpha,P);
		element_mul(aki2[i],ti1[i],Y);
		element_add(aki1[i],aki1[i],aki2[i]);
		element_mul(aki2[i],ti1[i],P);
		element_mul(aki3[i],ti1[i],Q1);
		element_mul(aki4[i],ti2[i],Q2);
		element_add(aki4[i],aki4[i],aki3[i]);

		element_mul(aki3[i],ti2[i],P);
	}
	

	startpart = pbc_get_time();
//4.2 User Report Generation
	//step-1

	for(int i=0;i<number_of_iot;i++){
		//mpz_init(di[i]);
		//mpz_init(ri[i]);
		//_mpz_random(di[i],512);
		//_mpz_random(ri[i],512);
		//mpz_init_set_ui(Ci[i],1);
		mpz_add_ui(di[i],di[i],i);
		mpz_add_ui(ri[i],ri[i],i);
		for(int j=0;j<number_of_type;j++){
			mpz_powm(temp_mpz,gi[j],di[i],n_square);
			mpz_mul(Ci[i],Ci[i],temp_mpz);
		}
		mpz_powm(temp_mpz,ri[i],n,n_square);
		mpz_mul(Ci[i],Ci[i],temp_mpz);
		mpz_mod(Ci[i],Ci[i],n_square);
		//gmp_printf("C%d : %Zd\n",i,Ci[i]);	
	}
	//step-2
	t sigmai[number_of_iot];
	for(int i=0;i<number_of_iot;i++){
		IG1(sigmai[i]);
		element_from_hash(sigmai[i],xi[i],100);
		element_mul(sigmai[i],xi[i],sigmai[i]);
		//element_printf("sigmai%d : %B\n",i,sigmai[i]);
	}
	endpart = pbc_get_time();
	user2 = pbc_get_time();
	user3 = user2 -user1;
//4.3 Privacy-Preserving Report Aggregation
	node1 = pbc_get_time();
	startpart=pbc_get_time();
	for(int i=0;i<number_of_iot;i++){
		element_add(temp_IG1,temp_IG1,sigmai[i]);

		element_from_hash(temp_IG1,sigmai[i],100);
		element_pairing(temp_IGT,Yi[i],temp_IG1);
		element_mul(temp_IGT,temp_IGT,temp_IGT);
	}
		element_pairing(temp_IGT,P,temp_IG1);
	//step-1
	mpz_t C;
	mpz_init_set_ui(C,1);
	for(int i=0;i<number_of_iot;i++){
		mpz_mul(C,C,Ci[i]);
		mpz_mod(C,C,n_square);
	}
	//step-2
	t sigmag;
	IG1(sigmag);
	element_from_hash(temp_IG1,C,100);
	element_mul(sigmag,xg,temp_IG1);
	endpart = pbc_get_time();
	agg3 = endpart-startpart;
	node2 = pbc_get_time();
	node3 = node3+node2 -node1;
// 4.4 Secure Report Reading and Response
	startpart=pbc_get_time();
	server1 = pbc_get_time();
	element_pairing(temp_IGT,P,sigmag);
	element_from_hash(temp_IG1,C,100);
	element_pairing(temp_IGT,Yg,sigmag);
	//step-1
	mpz_t M;
	mpz_init(M);
	mpz_set_ui(M,1);
	for(int i=0;i<number_of_type;i++){
		for(int j=0;j<number_of_type;j++){
			mpz_add(di[j],di[j],di[j]);
			mpz_mod(di[i],di[j],n);
		}
		mpz_mul(temp_mpz,ai[i],di[i]);
		mpz_mod(temp_mpz,temp_mpz,n);
		mpz_add(M,temp_mpz,M);
		mpz_mod(M,M,n);
//		gmp_printf("M%d:%Zd\n",i,M);
	}
	//step-2

		mpz_set(Xj[number_of_type-1],M);
	for(int j=number_of_type-1;j>=2;j--){
	//	mpz_init(Xj[j-1]);
	//	mpz_init(Dj[j]);
		mpz_mod(Xj[j-1],Xj[j],ai[j]);
		mpz_add_ui(Xj[j-1],Xj[j-1],j);
	//	gmp_printf("X%d:%Zd\n",j,Xj[j]);
		mpz_sub(Dj[j],Xj[j],Xj[j-1]);
		mpz_cdiv_q(Dj[j],Dj[j],ai[0]);
		//gmp_printf("D%d:%Zd\n",j,Dj[j]);
	}
	endpart = pbc_get_time();
	de3=endpart-startpart;
//Response
//step-1
	startpart = pbc_get_time();
	element_pow_zn(C1,ePPalpha,s);
	element_mul(C1,C1,m);

	element_mul(C2,s,P);
	
	element_mul(C3,s,Y);
	element_mul(C4,s,Q1);
	element_neg(C4,C4);
	element_sub(C3,C3,C4);

	element_mul(C4,s,Q2);
	element_neg(C4,C4);

	element_from_hash(temp_IG1,C1,100);
	element_mul(temp_IG1,temp_IG1,x);

	server2 = pbc_get_time();
	server3 = server3 + server2-server1;
	//step-2
	node1 = pbc_get_time();
	element_from_hash(temp_IG1,C1,100);
	element_pairing(temp_IGT,temp_IG1,Y);
	element_pairing(temp_IGT,P,sigma);
	node2 = pbc_get_time();
	node3 = node3+node2-node1;

	//Decryption of User
	user1 = pbc_get_time();
	for(int i=0;i<number_of_iot;i++){
		element_pairing(up,C2,aki1[i]);
		element_pairing(down,C3,aki2[i]);
		element_pairing(temp_IGT,aki3[i],C4);
		element_mul(down,temp_IGT,down);

		element_pairing(temp_IGT,aki4[i],C2);
		element_mul(down,temp_IGT,down);

		element_div(temp_IGT,up,down);
		element_div(temp_IGT,C1,temp_IGT);
	}
	user2 = pbc_get_time();
	user3 = user2 - user1 + user3;
	endpart = pbc_get_time();
	end = pbc_get_time();
	iot3 = user3;	
	printf("number:%d\n",number_of_iot);
	printf("total:%.3f\n",(end-start)*1000);
	printf("iot:%.3f\n",iot3*1000/number_of_iot);
	printf("node:%.3f\n",node3*1000);
	printf("server:%.3f\n",server3*1000);
	printf("agg:%.3f\n",agg3*1000);
	printf("de:%.3f\n",de3*1000);

}



