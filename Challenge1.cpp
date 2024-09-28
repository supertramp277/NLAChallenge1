#include <Eigen/Dense>
#include <iostream>
#include <random>
#include <unsupported/Eigen/SparseExtra>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using namespace Eigen;

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
    MatrixXd original_image_matrix(height, width);

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
    MatrixXd noised_image_matrix(height, width);
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

    // By map you can creat a vector reference to memeory without copying data
    Map<VectorXd> v(original_image_matrix.data(), original_image_matrix.size());
    Map<VectorXd> w(noised_image_matrix.data(), noised_image_matrix.size());

    // Verify the size of the vectors
    std::cout << "Original image vector v's size: " << v.size() << std::endl;
    std::cout << "Noisy image vector w's size: " << w.size() << std::endl;
    std::cout << "Euclidean norm of v is: " << v.norm() << std::endl;

    // Free memory
    stbi_image_free(image_data);

    return 0;
}