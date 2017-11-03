
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
Rcpp::List sortIdxByCol2(mat score, umat label, uword k) {
  //mat M(3, 2, fill::randu);
  //uvec indices = find(M > 0.5);
  uword nr = score.n_rows;
  uword nc = score.n_cols;
  
  //mat idxMat(nr, nc);
  //idxMat.each_col();
  
  umat idxMat2(nr, nc);
  
  mat sortedScore(nr, nc);
  
  umat sortedLabel(nr, nc);
  

  
  for (uword i = 0; i < nc; i++) {
    
    vec curCol = score.col(i);
    
    uvec colIdx = sort_index(curCol, "descend");
    
    
    // must '+ 1' for output into R
    //colIdx = colIdx + 1;
    //##Rcout << colIdx << "\n";
    //##vec colIdx_ = conv_to<vec>::from(colIdx);
    //##Rcout << colIdx_ << "\n";
    //##idxMat.col(i) = colIdx_;
    //idxMat2.col(i) = colIdx;
    
    sortedScore.col(i) = curCol(colIdx);
    
    uvec curL = label.col(i);
    
    sortedLabel.col(i) = curL(colIdx);
    
    
    
  }
  
  sortedScore = sortedScore.rows(0, k - 1);
  sortedLabel = sortedLabel.rows(0, k - 1);
  
  return Rcpp::List::create(Rcpp::Named("scoreSorted") = sortedScore,
                            Rcpp::Named("labelSorted") = sortedLabel);
  //umat t = ind2sub( size(M), indices );
  
  //Rcout << M << "\n";
  //Rcout << idxMat2 << "\n";
  
  //return idxMat2;
  //return sortedScore;
}



