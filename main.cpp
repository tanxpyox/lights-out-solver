#include <algorithm>
#include <cstddef>
#include <cassert>
#include <cmath>
#include <iostream>
#include <cstring>
 
#define H 5  //height of array
#define L 5 //length of array

template<typename MatrixType> struct matrix_traits {
  typedef typename MatrixType::index_type index_type;
  typedef typename MatrixType::value_type value_type;
  static index_type min_row(MatrixType const& A){ 
  	return A.min_row(); }
  static index_type max_row(MatrixType const& A){ 
  	return A.max_row(); }
  static index_type min_column(MatrixType const& A){ 
  	return A.min_column(); }
  static index_type max_column(MatrixType const& A){ 
  	return A.max_column(); }
  static value_type& element(MatrixType& A, index_type i, index_type k){ 
  	return A(i,k); }
  static value_type element(MatrixType const& A, index_type i, index_type k){ 
  	return A(i,k); }
};

// specialization of the matrix traits for built-in two-dimensional
// arrays
template<typename T, std::size_t rows, std::size_t columns>
 struct matrix_traits<T[rows][columns]> {
  typedef std::size_t index_type;
  typedef T value_type;
  static index_type min_row(T const (&)[rows][columns]) { 
  	return 0; }
  static index_type max_row(T const (&)[rows][columns]) { 
  	return rows-1; }
  static index_type min_column(T const (&)[rows][columns]) { 
  	return 0; }
  static index_type max_column(T const (&)[rows][columns]) { 
  	return columns-1; }
  static value_type& element(T (&A)[rows][columns],
                             index_type i, index_type k) { 
	return A[i][k]; }
  static value_type element(T const (&A)[rows][columns],
                            index_type i, index_type k) { 
	return A[i][k]; }
};

// Swap rows i and k of a matrix A
// Note that due to the reference, both dimensions are preserved for
// built-in arrays
template<typename MatrixType>
 void swap_rows(MatrixType& A,
                 typename matrix_traits<MatrixType>::index_type i,
                 typename matrix_traits<MatrixType>::index_type k) {
  matrix_traits<MatrixType> mt;
  typedef typename matrix_traits<MatrixType>::index_type index_type;

  // check indices
  assert(mt.min_row(A) <= i);
  assert(i <= mt.max_row(A));

  assert(mt.min_row(A) <= k);
  assert(k <= mt.max_row(A));

  for (index_type col = mt.min_column(A); col <= mt.max_column(A); ++col)
    std::swap(mt.element(A, i, col), mt.element(A, k, col));
}

// divide row i of matrix A by v
template<typename MatrixType>
 void divide_row(MatrixType& A,
                  typename matrix_traits<MatrixType>::index_type i,
                  typename matrix_traits<MatrixType>::value_type v) {
  matrix_traits<MatrixType> mt;
  typedef typename matrix_traits<MatrixType>::index_type index_type;

  assert(mt.min_row(A) <= i);
  assert(i <= mt.max_row(A));

  assert(v != 0);

  for (index_type col = mt.min_column(A); col <= mt.max_column(A); ++col)
    mt.element(A, i, col) /= v;
}

// in matrix A, add v times row k to row i
template<typename MatrixType>
 void add_multiple_row(MatrixType& A,
                  typename matrix_traits<MatrixType>::index_type i,
                  typename matrix_traits<MatrixType>::index_type k,
                  typename matrix_traits<MatrixType>::value_type v) {
  matrix_traits<MatrixType> mt;
  typedef typename matrix_traits<MatrixType>::index_type index_type;

  assert(mt.min_row(A) <= i);
  assert(i <= mt.max_row(A));

  assert(mt.min_row(A) <= k);
  assert(k <= mt.max_row(A));

  for (index_type col = mt.min_column(A); col <= mt.max_column(A); ++col)
    mt.element(A, i, col) = (abs(mt.element(A, i, col)) + v * abs(mt.element(A, k, col)))%2;
}

// convert A to reduced row echelon form
template<typename MatrixType>
 void to_reduced_row_echelon_form(MatrixType& A) {
  matrix_traits<MatrixType> mt;
  typedef typename matrix_traits<MatrixType>::index_type index_type;

  index_type lead = mt.min_row(A);

  for (index_type row = mt.min_row(A); row <= mt.max_row(A); ++row) {
    if (lead > mt.max_column(A))
      return;
    index_type i = row;
    while (mt.element(A, i, lead) == 0)     {
      ++i;
      if (i > mt.max_row(A)){
        i = row;
        ++lead;
        if (lead > mt.max_column(A))
          return;
      }
    }
    swap_rows(A, i, row);
    divide_row(A, row, mt.element(A, row, lead));
    for (i = mt.min_row(A); i <= mt.max_row(A); ++i) {
      if (i != row)
        add_multiple_row(A, i, row, mt.element(A, i, lead));
    }
  }
}

int main() {
  int M[H*L][H*L+1];
  std::memset(M,0,sizeof(M));
   
  for (int i=0;i<H;i++){
  	for (int j=0;j<L;j++){
  		int numbernow = i*L+j;
  		M[numbernow][numbernow] =1;
  		if (numbernow>=L) M[numbernow-L][numbernow] =1;        //if not on top top, toggle cell above
  		if (numbernow%L!=0) M[numbernow-1][numbernow] =1;      //if not on left edge, toggle left cell
  		if (numbernow%L!=L-1) M[numbernow+1][numbernow] =1;    //if not on right edge, toggle right cell
  		if (numbernow+L<=H*L-1) M[numbernow+L][numbernow] = 1; //if not on bottom edge, toggle cell below
	  }
  }
  
//input question matrix.
//0 =  target colour (to make the whole board the colour of 0); 1 = colour to eliminate  
  std::cout << "Enter a " << H << "*" << L << " matrix: (1 = colour to eliminate, 0 = keep)" << std::endl;
  for (int i=0;i<H;i++){
  	for (int j=0;j<L;j++){
  		std::cin >> M[i*L+j][H*L];
	  }
  }

  
  to_reduced_row_echelon_form(M);
/*
  int check = 0;
    for (int j = 0; j < H*L+1; ++j)
      check+= abs(M[H*L-1][j]); 
  if (check ==0) {
  	std::cout << "No solution!!";
  	return 0;
  }*/

//output
  std::cout << "Answer: (Please click on the cells marked with *)" << std::endl;
  for (int i = 0; i < H; ++i){
    for (int j = 0; j < L; ++j){
    	char output;
      	if (abs(M[i*L+j][H*L])) output = '*'; 
	else output = '.';
        std::cout << output  << ' ';
    }
    std::cout << "\n";
  }
 
  return EXIT_SUCCESS;
}
