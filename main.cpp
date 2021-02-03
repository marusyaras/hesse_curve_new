#include <iostream>
#include "tommath.h"
mp_err err;
struct coord{
    mp_int x;
    mp_int y;
    mp_int z;
};
void put_in_new_dot_val (coord* c, mp_int* x, mp_int* y, mp_int* z)
{
    err = mp_init (&(c->x));
    err = mp_init (&(c->y));
    err = mp_init (&(c->z));
    err = mp_copy (x, &(c->x));
    err = mp_copy (y, &(c->y));
    err = mp_copy (z, &(c->z));
}
void put_in_dot_val (coord* c, mp_int* x, mp_int* y, mp_int* z)
{
    err = mp_copy (x, &(c->x));
    err = mp_copy (y, &(c->y));
    err = mp_copy (z, &(c->z));
}

void addition(coord c1, coord c2, coord* c3, mp_int modulo){
    mp_int t1, t2, t3, t4, t5, t6, t7;
    err = mp_init_copy (&t1, &(c1.x));
    err = mp_init_copy (&t2, &(c1.y));
    err = mp_init_copy (&t3, &(c1.z));
    err = mp_init_copy (&t4, &(c2.x));
    err = mp_init_copy (&t5, &(c2.y));
    err = mp_init_copy (&t6, &(c2.z));
    err = mp_init (&t7);
    //T1 = X1
    //T2 = Y1
    //T3 = Z1
    //T4 = X2
    //T5 = Y2
    //T6 = Z2
    //T7 = T1*T6
    err = mp_mul (&t1, &t6, &t7);
    //T1 = T1*T5
    err = mp_mul (&t1, &t5, &t1);
    //T5 = T3*T5
    err = mp_mul (&t3, &t5, &t5);
    //T3 = T3*T4
    err = mp_mul (&t3, &t4, &t3);
    //T4 = T2*T4
    err = mp_mul (&t2, &t4, &t4);
    //T2 = T2*T6
    err = mp_mul (&t2, &t6, &t2);
    //T6 = T2*T7
    err = mp_mul (&t2, &t7, &t6);
    //T2 = T2*T4
    err = mp_mul (&t2, &t4, &t2);
    //T4 = T3*T4
    err = mp_mul (&t3, &t4, &t4);
    //T3 = T3*T5
    err = mp_mul (&t3, &t5, &t3);
    //T5 = T1*T5
    err = mp_mul (&t1, &t5, &t5);
    //T1 = T1*T7
    err = mp_mul (&t1, &t7, &t1);
    //T1 = T1-T4
    err = mp_sub (&t1, &t4, &t1);

    //T2 = T2-T5
    err = mp_sub (&t2, &t5, &t2);
    //T3 = T3-T6
    err = mp_sub (&t3, &t6, &t3);


    put_in_dot_val(c3, &t2, &t1, &t3);
    err = mp_mod(&(c3->x), &modulo, &(c3->x));
    err = mp_mod(&(c3->y), &modulo, &(c3->y));
    err = mp_mod(&(c3->z), &modulo, &(c3->z));


    mp_clear (&t1);
    mp_clear (&t2);
    mp_clear (&t3);
    mp_clear (&t4);
    mp_clear (&t5);
    mp_clear (&t6);
    mp_clear (&t7);


}

void montgomery_ladder (coord c, coord* res, mp_int modulo, mp_int numb){

    coord new_c;

    int am_bits = mp_count_bits (&numb);

    struct coord P, Q;
    mp_int temp, zero;
    err = mp_init_multi(&temp, &zero, NULL);
    mp_set (&zero, 0);

    mp_int x, y, z;
    err = mp_init_multi(&x, &y, &z, NULL);
    mp_set (&x, 1);
    err = mp_read_radix(&y, "-1", 10);
    mp_set (&z, 0);
    put_in_new_dot_val(&Q, &x, &y, &z);
    put_in_new_dot_val(&P, &(c.x), &(c.y), &(c.z));
    for (int i = am_bits-1; i>0; i--){
        mp_set (&temp, 1<<i);
        err = mp_and(&numb, &temp, &temp);
        if (mp_cmp(&temp, &zero)!=MP_EQ){
            addition(Q, P, &Q, modulo);
            addition(P, P, &P, modulo);
        }
        else{
            addition(P, Q, &P, modulo);
            addition(Q, Q, &Q, modulo);
        }
    }
    put_in_dot_val(res, &(Q.x), &(Q.y), &(Q.z));
    mp_clear(&x);
    mp_clear(&y);
    mp_clear(&z);
    mp_clear(&temp);
    mp_clear(&zero);
}

bool check_on_curve (coord* c, mp_int* modulo){
    mp_int x3, y3, z3, d, eq1, eq2, var;
    err = mp_init_multi (&x3, &y3, &z3, &d, &eq1, &eq2, &var, NULL);
    mp_set(&d, 3);
    err = mp_sqr(&(c->x), &x3);
    err = mp_mul(&x3, &(c->x), &x3);//xˆ3

    err = mp_sqr(&(c->y), &y3);
    err = mp_mul(&y3, &(c->y), &y3);//yˆ3

    err = mp_sqr(&(c->z), &z3);
    err = mp_mul(&z3, &(c->z), &z3);//zˆ3

    err = mp_add(&x3, &y3, &var);
    err = mp_add(&var, &z3, &eq1);//xˆ3 + yˆ3 + zˆ3

    err = mp_mul (&(c->x), &(c->y), &var);
    err = mp_mul (&var, &(c->z), &var);
    err = mp_mul(&var, &d, &var);
    mp_int koef;
    err = mp_init (&koef);
    mp_set (&koef, 3);
    err = mp_mul (&var,  &koef, &eq2); //3*d*x*y*z

    err = mp_mod(&eq1, modulo, &eq1);
    err = mp_mod(&eq2, modulo, &eq2);

    if (mp_cmp (&eq1, &eq2) == MP_EQ){return 1;}
    else return 0;


}
bool equal_dots (coord* c1, coord* c2, mp_int modulo){
    mp_int x1, x2, y1, y2, z1, z2;
    err = mp_init_multi(&x1, &x2, &y1, &y2, &z1, &z2, NULL);
    err = mp_copy (&(c1->x), &x1);
    err = mp_copy (&(c2->x), &x2);
    err = mp_copy (&(c1->y), &y1);
    err = mp_copy (&(c2->y), &y2);
    err = mp_copy (&(c1->z), &z1);
    err = mp_copy (&(c2->z), &z2);
    err = mp_mod (&x1, &modulo, &x1);
    err = mp_mod (&x2, &modulo, &x2);
    err = mp_mod (&y1, &modulo, &y1);
    err = mp_mod (&y2, &modulo, &y2);
    err = mp_mod (&z1, &modulo, &z1);
    err = mp_mod (&z2, &modulo, &z2);
    bool eqq = ((mp_cmp(&x1, &x2)==MP_EQ)&&(mp_cmp(&y1, &y2)==MP_EQ)&&(mp_cmp(&z1, &z2)==MP_EQ));
    mp_clear(&x1);
    mp_clear(&x2);
    mp_clear(&y1);
    mp_clear(&y2);
    mp_clear(&z1);
    mp_clear(&z2);
    return eqq;
}

int main() {
    mp_int p, x, y, z, d, m;
    err = mp_init(&p);
    err = mp_init(&d);
    err = mp_init(&x);
    err = mp_init(&y);
    err = mp_init(&z);
    err = mp_init(&m);
    err = mp_read_radix(&p, "115792089237316195423570985008687907853269984665640564039457584007913111864739",
                  10);
    err = mp_read_radix(&x, "93528601524582384978654134272795468222405403055717890280271688132874849008326", 10);
    err = mp_read_radix(&y, "14443324612566128911211262381388707474030458136470034119105598903952521080679", 10);
    mp_set (&d, 3);
    mp_set(&z, 1);
    err = mp_read_radix(&m, "115792089237316195423570985008687907852907080286716537199505774922796921406320", 10);
    struct coord c;
    put_in_new_dot_val(&c, &x, &y, &z);

    struct coord zero, time;
    mp_int zero_x, zero_y, zero_z;
    err = mp_init_multi (&zero_x, &zero_y, &zero_z, NULL);
    mp_set(&zero_x, 1);
    err = mp_read_radix(&zero_y, "-1", 10);
    mp_set(&zero_z, 0);
    put_in_new_dot_val(&zero, &zero_x, &zero_y, &zero_z);
    put_in_new_dot_val(&time, &zero_x, &zero_y, &zero_z);

    int bits_numb_m = mp_count_bits(&m);
    mp_int rand_num;
    err = mp_init (&rand_num);
    err = mp_rand (&rand_num, bits_numb_m-1);

    std::cout<< "Test1\n";
    montgomery_ladder(c, &zero, p, rand_num);
    if (check_on_curve(&zero, &p)){
        std::cout<<"Dot lays on curve\n";
    }
    else { std::cout<<"Dot does not lay on curve\n";}


    struct coord c2;
    mp_set(&zero_x, 1);
    err = mp_read_radix(&zero_y, "-1", 10);
    mp_set(&zero_z, 0);
    put_in_dot_val(&zero, &zero_x, &zero_y, &zero_z);
    put_in_new_dot_val(&c2, &zero_x, &zero_y, &zero_z);
    std::cout<<"Test2\n";
    montgomery_ladder(c, &c2, p, m);
    if (mp_cmp(&(c2.z), &zero_z)==MP_EQ){
        std::cout<<"OK\n";
    }
    else {std::cout<<"Not OK\n";}

    std::cout<<"Test3\n";
    err = mp_copy (&m, &rand_num);
    err = mp_add (&rand_num, &zero_x, &rand_num);//rand_num = m+1
    montgomery_ladder(c, &c2, p, rand_num);

    mp_int test1, test2;
    err = mp_init_multi(&test1, &test2, NULL);

    err = mp_mul (&(c2.x), &(c.y), &test1);
    err = mp_mul (&(c2.y), &(c.x), &test2);
    if (mp_cmp(&test1, &test2)==MP_EQ){
        std::cout<<"OK\n";
    }
    else {std::cout<<"Not OK\n";}

    err = mp_sub (&m, &zero_x, &rand_num);
    montgomery_ladder(c, &c2, p, rand_num);
    err = mp_mul (&(c2.x), &(c.x), &test1);
    err = mp_mul (&(c2.y), &(c.y), &test2);
    if (mp_cmp(&test1, &test2)==MP_EQ){
        std::cout<<"OK\n";
    }
    else {std::cout<<"Not OK\n";}

    std::cout<<"Test4\n";
    mp_int k, k1, k2;
    err = mp_init_multi(&k1, &k2, &k, NULL);
    err = mp_rand (&k1, 5);
    err = mp_rand (&k2, 5);
    err = mp_add (&k1, &k2, &k);
    struct coord r1, r2, r;

    put_in_new_dot_val(&r, &k, &k, &k);
    put_in_new_dot_val(&r1, &k, &k, &k);
    put_in_new_dot_val(&r2, &k, &k, &k);

    montgomery_ladder(c, &r1, p, k1);
    montgomery_ladder(c, &r2, p, k2);
    montgomery_ladder(c, &r, p, k);

    addition (r1, r2, &r1, p);
    err = mp_mul (&(r.x), &(r1.y), &test1);
    err = mp_mul (&(r.y), &(r1.x), &test2);

    if (mp_cmp(&test1, &test2)==MP_EQ){
        std::cout<<"OK\n";
    }
    else {std::cout<<"Not OK\n";}
    mp_clear(&k);
    mp_clear(&k1);
    mp_clear(&k2);
    mp_clear(&zero_x);
    mp_clear(&zero_y);
    mp_clear(&zero_z);
    mp_clear(&p);
    mp_clear(&z);
    mp_clear(&x);
    mp_clear(&y);
    mp_clear(&d);
    mp_clear(&m);


    return 0;

}
