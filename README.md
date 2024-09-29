# NLA Challenge1

- Use 256px image as input, print related matrix munually for checking, and output the genrated image.

- all the matrices are by rowmajor order and normalized to $[0,1]$ when processing. And convert to `Matrix<unsigned char>` when output.

- v.mtx represents v; w.mtx represents w. A1.mtx, A2.mtx, A3.mtx are repectively the $H_{av2}$, $H_{sh2}$ and $H_{lap}$ related convolution matrix.

- As we can see, there are some values which are greater than 1 or less than 0 for A2. So we can clip them to $[0,1]$ at the step of converting to `Matrix<unsigned char>` type.

- To ensure A's size is mn\*mn, we use zero padding here. And we store the generated sparse matrix by using triplet.
  ![Example Of Zero Padding](ZeroPadding.png)
