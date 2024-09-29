#include <Eigen/Dense>
#include <iostream>
#include <random>
#include <unsupported/Eigen/SparseExtra>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using namespace Eigen;
typedef Eigen::Triplet<double> T;

// Definite the convolution function by matrix production
SparseMatrix<double, RowMajor> convolution(const Matrix<double, Dynamic, Dynamic, RowMajor> &kernel, int height, int width)
{
    const int kernel_size = kernel.rows();
    const int m = height;
    const int n = width;
    const int mn = m * n;
    SparseMatrix<double, RowMajor> A(mn, mn);

    std::vector<T> hav2TripletList;
    hav2TripletList.reserve(mn * kernel_size * kernel_size);

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            int index_Ai = i * n + j; // row number of A1
            // ki, kj are respectively row index and column index of kernel, central point is (0,0)
            for (int ki = -kernel_size / 2; ki <= kernel_size / 2; ki++)
            {
                for (int kj = -kernel_size / 2; kj <= kernel_size / 2; kj++)
                {
                    int ci = i + ki; // contribute to n(width) shift each time when there is a vertical moving for convolution or inside the kernel
                    int cj = j + kj; // contribute just 1 shift each when there is a horizontal moving for convilution or inside the kernel
                    if (ci >= 0 && ci < m && cj >= 0 && cj < n)
                    {
                        int index_Aj = ci * n + cj; // column number of A1
                        hav2TripletList.push_back(Triplet<double>(index_Ai, index_Aj, kernel(ki + kernel_size / 2, kj + kernel_size / 2)));
                    }
                }
            }
        }
    }
    A.setFromTriplets(hav2TripletList.begin(), hav2TripletList.end());
    return A;
}

// Define the smooth function
void smoothImage(SparseMatrix<double, RowMajor> A1, int height, int width, VectorXd w)
{
    VectorXd smooth_image_vector = A1 * w;
    Matrix<double, Dynamic, Dynamic, RowMajor> smooth_image_matrix(height, width);
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            smooth_image_matrix(i, j) = smooth_image_vector(i * width + j);
        }
    }

    Matrix<unsigned char, Dynamic, Dynamic, RowMajor> smooth_image_output = smooth_image_matrix.unaryExpr(
        [](double pixel)
        { return static_cast<unsigned char>(pixel * 255); });
    const std::string smooth_image_path = "smoothedImage.png";
    stbi_write_png(smooth_image_path.c_str(), width, height, 1, smooth_image_output.data(), width);
}

// Define the sharpening function
void sharpenImage(const SparseMatrix<double, RowMajor> A2, int height, int width, const VectorXd v)
{
    // Perform a matrix production for sharpening and convert it to matrix
    VectorXd sharpen_image_vector = A2 * v;
    Matrix<double, Dynamic, Dynamic, RowMajor> sharpen_image_matrix(height, width);
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            sharpen_image_matrix(i, j) = sharpen_image_vector(i * width + j);
        }
    }

    // Convert the sharpened image to grayscale and export it using stbi_write_png
    Matrix<unsigned char, Dynamic, Dynamic, RowMajor> sharpen_image_output = sharpen_image_matrix.unaryExpr(
        [](double pixel)
        { return static_cast<unsigned char>(std::max(0.0, std::min(255.0, pixel * 255))); }); // ensure range [0,255]
    const std::string sharpen_image_path = "sharpenedImage.png";
    stbi_write_png(sharpen_image_path.c_str(), width, height, 1, sharpen_image_output.data(), width);
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <image_path>" << std::endl;
        return 1;
    }

    const char *input_image_path = argv[1];

    // Load the image using stb_image
    int width, height, channels;
    // for greyscale images force to load only one channel
    unsigned char *image_data = stbi_load(input_image_path, &width, &height, &channels, 1);
    if (!image_data)
    {
        std::cerr << "Error: Could not load image " << input_image_path << std::endl;
        return 1;
    }

    std::cout << "Image loaded: " << width << "x" << height << " with " << channels << " channels." << std::endl;

    // Convert the image_data to MatrixXd form, each element value [0,1]
    // Attention! if use MatrixXd here its columnmajor by default, we prefere rowmajor
    Matrix<double, Dynamic, Dynamic, RowMajor> original_image_matrix(height, width);

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            int index = (i * width + j) * channels;
            original_image_matrix(i, j) = static_cast<double>(image_data[index]) / 255;
        }
    }

    // Report the size of the matrix
    std::cout << "The Image Matrix Size Is: " << original_image_matrix.rows() << "*" << original_image_matrix.cols()
              << "=" << original_image_matrix.size() << std::endl;

    // Introduce noise and export
    Matrix<double, Dynamic, Dynamic, RowMajor> noised_image_matrix(height, width);
    // Uniform random number generator for noise
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(-50, 50);

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            int noise = distribution(generator);
            double noisedData = original_image_matrix(i, j) + static_cast<double>(noise) / 255;
            noised_image_matrix(i, j) = std::max(0.0, std::min(1.0, noisedData)); // ensure the value in [0,1]
        }
    }

    // Convert noise image to grayscale and export it by using stbi_write_png
    Matrix<unsigned char, Dynamic, Dynamic, RowMajor> noised_image_output = noised_image_matrix.unaryExpr(
        [](double pixel)
        { return static_cast<unsigned char>(pixel * 255); });
    const std::string noised_image_path = "NoisedImage.png";
    stbi_write_png(noised_image_path.c_str(), width, height, 1, noised_image_output.data(), width);

    // By map creat a vector reference to memeory without copying data, columnmajor by default,
    // however, we've already declared our data rowmajor so here rowmajor as well.
    Map<VectorXd> v(original_image_matrix.data(), original_image_matrix.size());
    Map<VectorXd> w(noised_image_matrix.data(), noised_image_matrix.size());

    // Verify the size of the vectors
    std::cout << "Original image vector v's size: " << v.size() << std::endl;
    std::cout << "Noisy image vector w's size: " << w.size() << std::endl;
    std::cout << "Euclidean norm of v is: " << v.norm() << std::endl;

    const int kernel_size = 3;
    // Create the kernel H_{av2}
    const double hav2_value = 1.0 / (kernel_size * kernel_size);
    Matrix<double, kernel_size, kernel_size, RowMajor> hav2;
    hav2.setConstant(hav2_value);

    // Perform convolution function
    SparseMatrix<double, RowMajor> A1 = convolution(hav2, height, width);
    // Check the nonzero numbers
    std::cout << "A1 nonzero numbers is " << A1.nonZeros() << std::endl;
    // Smooth the noise image
    smoothImage(A1, height, width, w);

    // Create the kernel H_{sh2}}
    Matrix<double, kernel_size, kernel_size, RowMajor> hsh2;
    hsh2 << 0.0, -3.0, 0.0,
        -1.0, 9.0, -3.0,
        0.0, -1.0, 0.0;
    // Perform convolution function
    SparseMatrix<double, RowMajor> A2 = convolution(hsh2, height, width);
    // Check the nonzero numbers
    std::cout << "A2 nonzero numbers is " << A2.nonZeros() << std::endl;
    // Verify if the matrix is symmetric, attention that A2.tranpose() should be type declaration
    double norm_diff = (A2 - SparseMatrix<double, RowMajor>(A2.transpose())).norm();
    std::cout << "A2 row:" << A2.rows() << "columns:" << A2.cols() << std::endl;
    std::cout << "Check if A2 is symmectric by its difference with transpose by norm:"
              << norm_diff << std::endl;
    // Sharpen the original image
    sharpenImage(A2, height, width, v);

    // Free memory
    stbi_image_free(image_data);

    return 0;
}