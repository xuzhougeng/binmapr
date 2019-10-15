#include <Rcpp.h>
using namespace Rcpp;

// Find window

// [[Rcpp::export]]
NumericMatrix findWindow(NumericVector pos, int windowLen) {

  int last_index = pos.length() - 1;

  // The vector only have one window

  if (pos[last_index] - pos[0] < windowLen){
    NumericMatrix index_mt(1,2) ;
    index_mt(0,0) = 1 ;
    index_mt(0,1) = last_index + 1;
    return index_mt;
  }

  //index matrix stores the start and end of each window
  int index =  0 ;
  int start_index = 0;
  int end_index = 0 ;
  NumericMatrix index_mt(last_index, 2);

  //phyiscal position of genotype
  int start_pos = pos[0] ;
  int end_pos = 0 ;

  // While loop
  int ptr = 0;

  while(ptr <= last_index){

    end_pos = pos[ptr] ;
    if ((end_pos - start_pos) > windowLen){

      //If middle position is isolated, the index of start
      //and end is same
      if (end_index < start_index){
        end_index = start_index ;
      }

      //fill in the matrix with the index of start and end
      index_mt(index, 0) = start_index + 1 ;
      index_mt(index, 1) = end_index + 1;
      //increase the matriw row index
      index ++ ;

      //set the start position
      start_index = ptr ;
      start_pos = end_pos ;

      // increase the ptr
      ptr ++ ;

    } else{
      end_index = ptr ;
      ptr = ptr + 1 ;

    }
  }

  //For the case that the last position is still in the window
  if (pos[last_index] - pos[start_index] <=  windowLen){
    index_mt(index,0) = start_index + 1;
    index_mt(index,1) = last_index  + 1;
  }

  // For the case that the last position is isolated
  if (pos[last_index] - pos[last_index - 1] > windowLen){
    index_mt(index, 0) = last_index + 1;
    index_mt(index, 1) = last_index + 1;
  }

  NumericMatrix mt = index_mt( Range(0, index), Range(0,1));

  return mt;

}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
geno <- c(1,10:15, 20, 30:38, 51,61)
findWindow(geno, windowLen = 10)
*/
