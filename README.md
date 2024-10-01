# NLA Challenge 1

Hands-on Challenge 1 of the course _Numerical Linear Algebra_ by Professor Antonietti, Polimi, a.y. 2024/25.

## Execution

- To run the main code `Challenge1.cpp` on terminal:

  ```bash
  g++ -I ${mkEigenInc} Challenge1.cpp -o exec
  ./exec einstein.jpg > output.txt
  ```

- This is the output file: [output.txt](output.txt)

## Explanation

- Use 256 pixel image _einstein.jpg_ as input, print related matrix manually for checking, and output the generated image.

- all the matrices are by rowmajor order and normalized to $[0,1]$ when processing. And convert to `Matrix<unsigned char>` when output.

- v.mtx represents $v$; w.mtx represents $w$. A1.mtx, A2.mtx, A3.mtx are repectively the $H_{av2}$, $H_{sh2}$ and $H_{lap}$ related convolution matrices.

- As we can see, there are some values which are greater than 1 or less than 0 for `A2`. So we can clip them to $[0,1]$ at the step of converting to `Matrix<unsigned char>` type.

- There are some common functions for convenience

  - `convolutionMatrix()`: Definite the convolution function by matrix product, return the sparse matrix of dimension $(mn,mn)$.
  - `convolutionMatrix2()`: Alternative function.
  - `void outputImage(const VectorXd &vectorData, int height, int width, const std::string &path)` Define the function that convert a **vector** to a `height`\*`width` rectangle matrix. Finally convert to `Matrix<unsigned char>` type and output it to _image.png_ by `stbi_write_png()` method.
  - `exportVector()`: export a vector to a mtx file. **Please pay attention here if we use `saveMarketVector()`, then when using LIS reading, some issues may occur**.
  - `exportSparsematrix():` export a sparse matrix to a mtx file by using `saveMarket()` method.
  - `isSymmetric():` check if a matrix is symmetric by comparing the transpose of it and see the norm of the difference of the two matrix. If the norm is less than the tolerance, we can say that the matrix is symmetric.
  - `isPositiveDefinite()`: check if a matrix is positive definite by `cholesky()` (searched online). So, if the matrix is symmetic positive definite, conjungate gradient solver can be used.

- To ensure $A$'s size is $(mn,mn)$, we use zero padding here. So, the pixel $A(0,0)=60$ is transformed into $A_{\text{new}}=60*5-73-113=114$.
  ![Example Of Zero Padding](ZeroPadding.png)
