#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

/* Decide the genotype in window with a simple
 rule.
 - NA Ratio > 0.75, genotype is consided as missing
 - The most frequent genotype is the represent
 - genotype 0 == genotype 1, 0
 - genotype 2 == genotype 1, 2
 - genotype 0 == genotype 2, 1
 - Otherwise: missing
 */
int decideGeno(int g0, int g1, int g2, int gna,
               int minSites, float na_ratio){

  float total = g0 + g1 + g2 + gna;

  if ( ((gna / total) > na_ratio)  & (total > minSites) ){
    return 3;
  }

  if ( (g0 > g1) & (g0 > g2) & (g0 > gna) ){
    return 0;
  } else if ( (g1 > g0) & (g1 > g2) & (g1 > gna)  ){
    return 1;
  } else if ( (g2 > g0) & (g2 > g1) & (g2 > gna)){
    return 2;
  } else if ( (gna > g0) & (gna > g1) & (gna > g2)) {
    return 3;
  } else if ( (g0 == g1) | (g0 == gna) ){
    return 0;
  } else if ( (g2 == g1) | (g2 == gna) ){
    return 2;
  } else if ( g0 == g2 ){
    return 1;
  } else {
    return 3;
  }
}


// [[Rcpp::export]]
NumericVector findRepresent(NumericMatrix mt, NumericVector geno){

  // initialize the output numeric vector
  int geno_num = mt.nrow();
  NumericVector new_x(geno_num);

  CharacterVector old_x_names = geno.names();
  CharacterVector new_x_names(geno_num);

  // initialize the variable in loop
  int window_start, window_end, i, j, temp;

  // initialize the genotype of 0, 1, 2 and NA
  unsigned geno_0 = 0, geno_1 = 0, geno_2 = 0, geno_NA = 0;

  // traverse the matrix and count the genotype in geno vector
  for (i = 0; i < geno_num; i++){

    window_start = mt(i,0);
    window_end = mt(i,1);

    // if the window only have one site
    if (window_start == window_end){
      new_x[i] = geno[window_start-1];
      new_x_names[i] = old_x_names[window_start-1];
      continue;
    }


    // Otherwise restore counter
    geno_0 = 0, geno_1 = 0, geno_2 = 0, geno_NA = 0;

    for (j = window_start - 1 ; j < window_end; j++){
      //count the NA, 1,0,2
      if (NumericVector::is_na(geno[j])){
        ++geno_NA;
        continue;
      }
      temp = int(geno[j]);
      switch(temp){
      case 0:
        ++geno_0;
        break;
      case 1:
        ++geno_1;
        break;
      case 2:
        ++geno_2;
        break;
      }
    }

    temp = decideGeno(geno_0, geno_1, geno_2, geno_NA, 10, 0.5);
    if (temp == 3){
      new_x[i] = NA_REAL;
      new_x_names[i] = old_x_names[window_start-1];
    } else{
      new_x[i] = temp;
      new_x_names[i] = old_x_names[window_start-1];
    }

  }
  new_x.names() = new_x_names;
  return new_x;
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
geno <- c(rep(1, 10),0, rep(1,3),2,rep(2,3))
names(geno) <- c(1,10:15, 20, 30:38, 51)
index_mt <- matrix(c(1,3,4,8,9,17,18,18), byrow = T, ncol = 2)
findRepresent(index_mt, geno)
*/
