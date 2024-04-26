#include <iostream>
#include <vector>



class Matrix {
  size_t cols;
  size_t rows;
  std::vector<std::vector<double>> data;
  public:
  Matrix(size_t rows, size_t cols) : rows(rows), cols(cols) {
    data.resize(rows, std::vector<double>(cols));
  }
  double&  operator()(size_t i, size_t j) {return data[i][j];}
  void print() {
    size_t i, j;
    for(i = 0; i < rows; i++) {
      for(j = 0; j < cols; j++) {
        std::cout << data[i][j] << " ";
      }
    std::cout << std::endl;
    }
  }
};

int main() {

   Matrix A(2,2);
   A(0,0) = 1, A(1,1)=1;
   A.print();
  return 0;
}
