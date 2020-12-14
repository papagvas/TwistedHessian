#include <tommath.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

struct PointProjective {
  mp_int* x;
  mp_int* y;
  mp_int* z;
};

int myNeg(mp_int* a, mp_int* p, mp_int* res) {
  mp_neg(a, a);
  mp_mod(a, p, res); 
  return 0;	
}

void mp_print(mp_int* outputpls, char* fancytext) {
  char buf[4096];                                      
  size_t length;                                                                                                     
  mp_to_radix(outputpls, buf, sizeof(buf), &length, 10);
  printf("%s", fancytext);  
  printf(" == %s, length = %zu \n", buf, length);  
}

void printPoint(struct PointProjective* whatever, char* text) {
  printf("%s\n", text);
  mp_print(whatever->x, "x");
  mp_print(whatever->y, "y");
  mp_print(whatever->z, "z");
}

void printAffine(struct PointProjective* proj, mp_int* p) {
  mp_int afx, afy, invz;
  mp_init_multi(&afx, &afy, &invz, NULL);

  mp_invmod(proj->z, p, &invz);
  mp_mulmod(proj->x, &invz, p, &afx);
  mp_mulmod(proj->y, &invz, p, &afy);

  mp_print(&afx, "affine x");
  mp_print(&afy, "affine y");
}

//void initPoint();
int exclZeroes(unsigned char* buf, mp_int* k) {
  int flag = 0;
  int ind = mp_ubin_size(k) * 8 - 1;
  int res = 0;
  while (flag == 0) {
    if (mp_get_bit(k, ind) == MP_YES)
      flag = 1;
    else {
      res = res + 1;
      ind =  ind - 1;
    }
      
  }
  return ind;
}

void convertToHesse(const mp_int* p, const mp_int* a, const mp_int* b, mp_int* res) {
  mp_int temp1, temp2, eight, four, x1, x2;
  const mp_int temp3, temp4, temp5, temp6, temp7, ten, finres; 
  mp_init_multi(&temp1, &temp2, &temp3, &temp4, &temp5, &temp6, &temp7, &eight, &four, &ten, &x1, &x2, &finres, NULL); //init temporary variables
  
  mp_read_radix(&temp2, "54", 10); // temp2 = 54
  mp_read_radix(&eight, "8", 10); // eight = 8
  mp_read_radix(&temp1, "400" ,10);
  mp_read_radix(&four, "4" , 10); 
  mp_read_radix(&ten, "10", 10);
  mp_read_radix(&temp7, "2", 10);
  mp_invmod(&temp2, p, &temp2); // temp3 = 54^(-1) (p)
  mp_exch(&temp2, &temp3);
  mp_mulmod(&temp3, b, p, &temp2); // temp2 = b * temp2 (p)
  mp_exch(&temp2, &temp4);
  mp_addmod(&temp4, &eight, p, &temp2); //temp2 = temp2 + 8 (p)
  mp_exch(&temp2, &temp5);
  myNeg(&temp2, p, &temp2); // temp2 = -temp2 (p)   
  mp_mulmod(&temp5, &four, p, &temp2); // 4*c
  mp_exch(&temp2, &temp6); 
  mp_addmod(&temp6, &temp1, p, &temp2); //temp1 = discr
  mp_sqrtmod_prime(&temp2, p, res); //res = sqrt(discr)
  myNeg(res, p, res);
  mp_exch(res, &finres);
  mp_invmod(&temp7, p, &temp7);
  mp_mulmod(&finres, &temp7, p, res); 
  mp_addmod(&ten, res, p, res); // D^3
  mp_exch(res, &finres);
  mp_read_radix(&temp7, "27", 10);  //105 = D(D^3 + 8)
  mp_invmod(&temp7, p, &temp7);
  mp_mulmod(a, &temp7, p, &temp2); // a / 27
  mp_addmod(&finres, &eight, p, &temp1); //D^3 + 8
  mp_invmod(&temp1, p, &temp1);
  mp_mulmod(&temp1, &temp2, p, res); //res = D
}

void nu(struct PointProjective* weierPoint, mp_int* d, mp_int* p, mp_int* res) {
  mp_int six, tree, one, nine, d_cubed, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, ts;
  mp_init_multi(&six, &d_cubed, &one, &nine, &ts, &temp1, &temp2, &tree, &temp3, &temp4, &temp5, &temp6, &temp7, &temp8, NULL);
  mp_read_radix(&six, "6", 10); //six 
  mp_read_radix(&one, "1", 10); // one
  mp_read_radix(&tree, "3", 10); // three
  mp_exptmod(d, &tree, p, &d_cubed); // D^3 
  mp_submod(&d_cubed, &one, p, &temp1); //D^3 - 1 = temp1
  mp_mulmod(&six, &temp1, p, &temp2); // temp2 = 6*(D^3 - 1)
  mp_read_radix(&nine, "9", 10); //nine
  mp_mulmod(&nine, &d_cubed, p, &temp3); // temp3 = 9D^3
  mp_addmod(weierPoint->y, &temp3, p, &temp4); // temp4 = y + 9D^3
  mp_mulmod(&tree, d, p, &temp5); // temp5 = 3*D
  mp_mulmod(&temp5, weierPoint->x, p, &temp6); // temp6 =  3*D*x
  mp_read_radix(&ts, "36", 10);
  mp_addmod(&temp6, &ts, p, &temp7); //temp7 = 3dx + 36
  mp_submod(&temp4, &temp7, p, &temp8); //temp8 = y + 9d^3 - 3dx - 36
  mp_mulmod(&temp8, &temp2, p, &temp8); // numerator of nu

  mp_mulmod(&d_cubed, &tree, p, &temp1); // temp1 = d^3 * 3
  mp_mulmod(d, weierPoint->x, p, &temp2); // temp2 = dx 
  mp_read_radix(&six, "12", 10); // six = 12 
  mp_addmod(&six, &temp2, p, &temp2); //dx + 12 
  mp_submod(&temp1, &temp2, p, &temp1); //3d^3 - dx - 12 
  mp_exptmod(&temp1, &tree, p, &temp1); // (3D^3 - Dx - 12)^3
  mp_read_radix(&six, "9", 10);
  mp_read_radix(&temp3, "2", 10);
  mp_exptmod(d, &temp3, p, &temp3); // temp3 = D^2 
  mp_mulmod(&six, &temp3, p, &temp3); // temp3 = temp3*9
  mp_addmod(weierPoint->x, &temp3, p, &temp3); //temp3 + x
  mp_exptmod(&temp3, &tree, p, &temp3); // temp3^3
  mp_addmod(&temp3, &temp1, p, &temp2); //temp3 + temp1
  mp_invmod(&temp2, p, &temp2); 
  mp_mulmod(&temp8, &temp2, p, res);
}

int isOnCurve(struct PointProjective* hessePoint, mp_int* a, mp_int* d, mp_int* p) {
  mp_int temp1, temp2, temp3, temp4, three, zero;
  mp_init_multi(&temp1, &temp2, &temp3, &temp4, &three, &zero, NULL);
  mp_read_radix(&zero, "0", 10);
  mp_read_radix(&three, "3", 10);
  mp_exptmod(hessePoint->x, &three, p, &temp1);
  mp_mulmod(&temp1, a, p, &temp1); // aX^3
  mp_exptmod(hessePoint->y, &three, p, &temp3);
  mp_addmod(&temp1, &temp3, p, &temp1); // aX^3 + Y^3
  mp_exptmod(hessePoint->z, &three, p, &temp4);
  mp_addmod(&temp1, &temp4 , p, &temp1); // aX^3 + Y^3 + Z^3
  mp_mulmod(hessePoint->x, hessePoint->y, p, &temp2); // x*y
  mp_mulmod(hessePoint->z, &temp2, p, &temp2); // x*y*z
  mp_mulmod(d, &temp2, p, &temp2); // d * x * y * z
  mp_submod(&temp1, &temp2, p, &temp1); 
  return mp_cmp_mag(&temp1, &zero) == MP_EQ;
}

void transformPoint(struct PointProjective* weierPoint, mp_int* nu, mp_int* d, mp_int* p, struct PointProjective* hesse) {
  mp_int nine, three, d_cubed, one, two, twelve, temp1, temp2;
  mp_init_multi(&nine, &three, &one, &two, &twelve, &d_cubed, &temp1, &temp2, NULL);
  mp_read_radix(&nine, "9", 10);
  mp_read_radix(&three, "3", 10);
  mp_read_radix(&one, "1", 10);
  mp_read_radix(&two, "2", 10);
  mp_read_radix(&twelve, "12", 10);
  mp_exptmod(d, &three, p, &d_cubed); // d_cubed = d^3
  mp_exptmod(d, &two, p, &temp1); //temp1 = d^2
  mp_mulmod(&nine, &temp1, p, &temp1); // temp1 = d^2*9
  mp_addmod(weierPoint->x, &temp1, p, &temp1); //temp1 = d^2 + 9 + x
  mp_mulmod(nu, &temp1, p, hesse->x); // x = u in hesse
  
  mp_mulmod(&three, &d_cubed, p, &temp1); // temp1 = 3d^3
  mp_mulmod(d, weierPoint->x, p, &temp2); // temp2 = d*x
  mp_addmod(&twelve, &temp2, p, &temp2); // temp1 = dx + 12
  mp_submod(&temp1, &temp2, p, &temp1); // temp1 = 3d^3 - dx - 12
  mp_mulmod(&temp1, nu, p, &temp2); //temp1 * nu 
  mp_submod(&temp2, &one, p, hesse->y); // y = v in hesse  

  mp_read_radix(hesse->z, "1", 10);
}

void standardSum(struct PointProjective* l, struct PointProjective* r, mp_int* p, struct PointProjective* res) {
  mp_int two, temp1, temp2, temp3, temp4;
  mp_init_multi(&two, &temp1, &temp2, &temp3, &temp4, NULL);
  mp_read_radix(&two, "2", 10);

  mp_exptmod(l->x, &two, p, &temp1);  // temp1 = x1^2
  mp_mulmod(r->y, r->z, p, &temp2); // temp2 = y2*z2 
  mp_mulmod(&temp1, &temp2, p, &temp1); // temp1 = temp1*temp2
  mp_exptmod(r->x, &two, p, &temp3); // temp3 = x2^2
  mp_mulmod(l->y, l->z, p, &temp4); // temp4 = y1*z1
  mp_mulmod(&temp3, &temp4, p, &temp3); // temp3 = temp3*temp4
  mp_submod(&temp1, &temp3, p, res->x); // res->x = temp1*temp3 

  mp_exptmod(l->z, &two, p, &temp1); // temp1 = 
  mp_mulmod(r->x, r->y, p, &temp2);
  mp_mulmod(&temp1, &temp2, p, &temp1);
  mp_exptmod(r->z, &two, p, &temp3);
  mp_mulmod(l->x, l->y, p, &temp4);
  mp_mulmod(&temp3, &temp4, p, &temp3);
  mp_submod(&temp1, &temp3, p, res->y);

  mp_exptmod(l->y, &two, p, &temp1);
  mp_mulmod(r->x, r->z, p, &temp2);
  mp_mulmod(&temp1, &temp2, p, &temp1);
  mp_exptmod(r->y, &two, p, &temp3);
  mp_mulmod(l->x, l->z, p, &temp4);
  mp_mulmod(&temp3, &temp4, p, &temp3);
  mp_submod(&temp1, &temp3, p, res->z);

}

void rotatedSum(struct PointProjective* l, struct PointProjective* r, mp_int* a, mp_int* p, struct PointProjective* res) { //rotated sum
  mp_int two, temp1, temp2, temp3, temp4;
  mp_init_multi(&two, &temp1, &temp2, &temp3, &temp4, NULL);
  mp_read_radix(&two, "2", 10);

  mp_exptmod(r->z, &two, p, &temp1);
  mp_mulmod(l->x, l->z, p, &temp2);
  mp_mulmod(&temp1, &temp2, p, &temp1);
  mp_exptmod(l->y, &two, p, &temp3);
  mp_mulmod(r->x, r->y, p, &temp4);
  mp_mulmod(&temp3, &temp4, p, &temp3);
  mp_submod(&temp1, &temp3, p, res->x);

  mp_exptmod(r->y, &two, p, &temp1);
  mp_mulmod(l->y, l->z, p, &temp2);
  mp_mulmod(&temp1, &temp2, p, &temp1);
  mp_exptmod(l->x, &two, p, &temp3);
  mp_mulmod(r->x, r->z, p, &temp4);
  mp_mulmod(&temp3, &temp4, p, &temp3);
  mp_mulmod(&temp3, a, p, &temp3);
  mp_submod(&temp1, &temp3, p, res->y);

  mp_exptmod(r->x, &two, p, &temp1);
  mp_mulmod(l->x, l->y, p, &temp2);
  mp_mulmod(&temp1, &temp2, p, &temp1);
  mp_mulmod(&temp1, a, p, &temp1);
  mp_exptmod(l->z, &two, p, &temp3);
  mp_mulmod(r->y, r->z, p, &temp4);
  mp_mulmod(&temp3, &temp4, p, &temp3);
  mp_submod(&temp1, &temp3, p, res->z);

}

void sum(struct PointProjective* l, struct PointProjective* r, mp_int* a, mp_int* p, struct PointProjective* res) {
  mp_int zero;
  mp_init(&zero);
  mp_read_radix(&zero, "0", 10);
  standardSum(l, r, p, res);
  if ((mp_cmp_mag(res->x, &zero) == MP_EQ) && (mp_cmp_mag(res->y, &zero) == MP_EQ) && (mp_cmp_mag(res->z, &zero) == MP_EQ))
    rotatedSum(l, r, a, p, res);
}

void montgomeryLadder(struct PointProjective* point, mp_int* k, mp_int* a, mp_int* p, struct PointProjective* res) {
  mp_int one, rx, ry, rz, tx1, ty1, tz1, tx2, ty2, tz2;
  mp_init_multi(&one, &rx, &ry, &rz, &tx1, &ty1, &tz1, &tx2, &ty2, &tz2, NULL);
  mp_read_radix(res->x, "0", 10);
  mp_read_radix(res->y, "1", 10);
  mp_read_radix(res->z, "1", 10);
  mp_read_radix(&one, "1", 10);
  myNeg(res->y, p, res->y);

  mp_int idx, idy, idz;
  mp_init_multi(&idx, &idy, &idz, NULL);
  
  mp_read_radix(&idx, "0", 10);
  mp_read_radix(&idy, "1", 10);
  mp_read_radix(&idz, "1", 10);
  myNeg(&idy, p, &idy);
  
  struct PointProjective id;
  id.x = &idx;
  id.y = &idy;
  id.z = &idz;

  struct PointProjective temp1;
  temp1.x = &tx1;
  temp1.y = &ty1;
  temp1.z = &tz1;

  struct PointProjective temp2;
  temp2.x = &tx2;
  temp2.y = &ty2;
  temp2.z = &tz2;


  mp_exch(point->x, &rx);
  mp_exch(point->y, &ry);
  mp_exch(point->z, &rz);
  
  struct PointProjective r;
  r.x = &rx;
  r.y = &ry;
  r.z = &rz;

 // printPoint(&r, "r");
//  printPoint(res, "q");

  unsigned char buf[mp_ubin_size(k)];
  mp_to_unsigned_bin(k, buf);// sizeof(buf), size);
  int ind = exclZeroes(buf, k);
  int iter = 1;
  for (int i = ind; i >= 0; i--) {
    if (mp_get_bit(k, i) == MP_NO) {
   //   printf("%d\n", iter);
    //  printf("%d\n", 0);
    //  sum(&r, res, a, p, &r);
      standardSum(&r, res, p, &temp1);
      rotatedSum(res, res, a, p, &temp2);
      sum(&temp1, &id, a, p, &r);
      sum(&temp2, &id, a, p, res);
    //  printAffine(&r, p);
    //  printAffine(res, p);
    }
    else {
    //  printf("%d\n", iter);
    //  printf("%d\n", 1);
    //  sum(res, &r, a, p, res);
      standardSum(res, &r, p, &temp1);
      rotatedSum(&r, &r, a, p, &temp2);
      sum(&temp1, &id, a, p, res);
      sum(&temp2, &id, a, p, &r);
 //     printAffine(&r, p);
  //    printAffine(res, p);

    }
    iter = iter + 1;
  //  printf("%d\n", i);
  }
  
}

int main() {
  mp_int three, p, neg_a, b, d, n, x, y, z, x1, y1, z1, m, test, d_hesse; //Weierstrass params
  mp_init_multi(&p, &neg_a, &b, &d, &n, &x, &y, &z, &x1, &y1, &z1, &m, &three, &test, &d_hesse, NULL); //Initialize them
  mp_read_radix(&p, "115792089237316195423570985008687907853269984665640564039457584007913111864739", 10); 
  mp_read_radix(&b, "9774",10);
  mp_read_radix(&neg_a, "2835",10); // -a
  
  convertToHesse(&p, &neg_a, &b, &d);
  mp_print(&d, "d");
  mp_read_radix(&x, "8119021769418921565894388378008753471190077835025535568042996557746564092404", 10);
  mp_read_radix(&y, "28771053522948763457620123068270166170894398578116419877580784537256379313025", 10);
  mp_read_radix(&z, "1", 10);
  
  struct PointProjective point;
  point.x = &x;
  point.y = &y;
  point.z = &z;
  nu(&point, &d, &p, &n);
  mp_print(&n, "nu"); 
  
  struct PointProjective hesse;
  hesse.x = &x1;
  hesse.y = &y1;
  hesse.z = &z1;
 
  transformPoint(&point, &n, &d, &p, &hesse); // transform to hessian form
  printPoint(&hesse, "hessian projective coordinates");
  printf("%s\n", "affine coordinates");
  printAffine(&hesse, &p);
 
  mp_read_radix(&three, "3", 10);
  mp_mulmod(&d, &three, &p, &d_hesse);
 
  mp_int idx, idy, idz, one, resx, resy, resz, test1, test2, test3, negx, negy, k1, k2, k3;
  mp_init_multi(&idx, &idy, &idz, &one, &resx, &resy, &resz, &test1, &test2, &test3, &negx, &negy, &k1, &k2, &k3, NULL);
  
  mp_read_radix(&one, "1", 10);
  mp_read_radix(&idx, "0", 10);
  mp_read_radix(&idy, "1", 10);
  mp_read_radix(&idz, "1", 10);
  myNeg(&idy, &p, &idy);
  
  if (isOnCurve(&hesse, &one, &d_hesse, &p))
    printf("%s\n", "Point is on the curve!");
  
  struct PointProjective id;
  id.x = &idx;
  id.y = &idy;
  id.z = &idz;

  struct PointProjective res;
  res.x = &resx;
  res.y = &resy;
  res.z = &resz;
 
  struct PointProjective testt;
  testt.x = &test1;
  testt.y = &test2;
  testt.z = &test3;

  transformPoint(&point, &n, &d, &p, &testt);
  mp_read_radix(&m, "115792089237316195423570985008687907852907080286716537199505774922796921406320", 10);
  printf("\n%s\n", "test2, mP = id");
  montgomeryLadder(&hesse, &m, &one, &p, &res);
  printAffine(&res, &p);
  
  printf("%s\n", "test3.1 [m+1]P = P");
  mp_addmod(&m, &one, &p, &m);
  montgomeryLadder(&testt, &m, &one, &p, &res);
  printAffine(&res, &p);

  transformPoint(&point, &n, &d, &p, &testt);
  printf("%s\n", "test3.2 [m-1]P = -P");
  printf("%s\n", "-p = (-x/y, 1/y)");
  mp_invmod(testt.y, &p, &negy);
  mp_mulmod(testt.x, &negy, &p, &negx);
  mp_print(&negx, "x of -P is");
  mp_print(&negy, "y of -P is");
  mp_submod(&m, &one, &p, &m);
  mp_submod(&m, &one, &p, &m); // m - 1
  montgomeryLadder(&testt, &m, &one, &p, &res);
  printAffine(&res, &p);
  printf("%s\n", "test 4 k1P + k2P = [k1+k2]P"); 
  mp_rand(&k1, 2);
  mp_rand(&k2, 2);
  mp_addmod(&k1, &k2, &p, &k3);
  transformPoint(&point, &n, &d, &p, &testt);
  transformPoint(&point, &n, &d, &p, &hesse);
  mp_addmod(&m, &one, &p, &m);
  mp_print(&k1, "k1");
  mp_print(&k2, "k2");
  montgomeryLadder(&hesse, &k1, &one, &p, &res); // res = k1P
  montgomeryLadder(&testt, &k2, &one, &p, &hesse); // hesse = k2P
  standardSum(&res, &hesse, &p, &id); // id = k1P + k2P
  printf("%s\n", "k1P + k2P");
  printAffine(&id, &p);
  transformPoint(&point, &n, &d, &p, &testt);
  montgomeryLadder(&testt, &k3, &one, &p, &res);
  printf("%s\n", "[k1+k2]P");
  printAffine(&res, &p);
  return 0;
}
