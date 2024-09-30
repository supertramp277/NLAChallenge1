# NLA Challenge1

- Use 256px image as input, print related matrix munually for checking, and output the genrated image.

- all the matrices are by rowmajor order and normalized to $[0,1]$ when processing. And convert to `Matrix<unsigned char>` when output.

- v.mtx represents v; w.mtx represents w. A1.mtx, A2.mtx, A3.mtx are repectively the $H_{av2}$, $H_{sh2}$ and $H_{lap}$ related convolution matrix.

- As we can see, there are some values which are greater than 1 or less than 0 for A2. So we can clip them to $[0,1]$ at the step of converting to `Matrix<unsigned char>` type.

- There are some common functions for convenience

  - `convolutionMatrix()`: Definite the convolution function by matrix production, return the sparsematrix mn\*mn.
  - `void outputImage(const VectorXd &vectorData, int height, int width, const std::string &path)` Define the function that convert a **vector** to a height\*width rectangle matrix. Finally convert to `Matrix<unsigned char>` type and output it to image.png by `stbi_write_png()` method.
  - `exportVector()`: export a vector to a mtx file.
  - `exprotSparsematrix():` export a sparse matrix to a mtx file by using `saveMarket()` method.
  - `isSymmetric():` check if a matrix is symmetric by comparing the transpose of it. And see the norm of the difference of the two matrix. If the norm is less than tolerance, we can say the matrix is symmetric.
  - `isPositiveDefinite()`: check if a matrix is positive definite by `cholesky()` (searched online). So if the matrix is symmetic positive definite, conjungate gradient solver can be used.

- To ensure A's size is mn\*mn, we use zero padding here. And we store the generated sparse matrix by using triplet.
  ![Example Of Zero Padding](ZeroPadding.png)
