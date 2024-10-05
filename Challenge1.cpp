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

SparseMatrix<double, RowMajor> convolution3(MatrixXd kernel, int height, int width,int kernel_size)
{      
    SparseMatrix<double, RowMajor> A(height*width, height*width);
    std::vector<T> hav2TripletList;
    hav2TripletList.reserve(height*width * kernel_size * kernel_size);
    //4 corner elements
    for (int m=1;m<kernel_size;m++ )
    {
      //top left
      hav2TripletList.push_back(Triplet<double>(0,0+m,kernel(1,m)));
      if(kernel(2,m)!=0){
      hav2TripletList.push_back(Triplet<double>(0,width+m-1,kernel(2,m)));}
      //top right
      hav2TripletList.push_back(Triplet<double>(width-1,width-1-m,kernel(1,m-1)));  
      if(kernel(2,m-1)!=0){
      hav2TripletList.push_back(Triplet<double>(width-1,2*width-2-m,kernel(2,m-1)));}
      //bottom left
      if(kernel(0,m)!=0){
      hav2TripletList.push_back(Triplet<double>(width*(height-1),width*(height-1)+m,kernel(0,m))); }
      hav2TripletList.push_back(Triplet<double>(width*(height-1),width*(height-1)+m-width-1,kernel(1,m)));  
      //bottom right
      if(kernel(0,m-1)!=0){
      hav2TripletList.push_back(Triplet<double>(width*height-1,width*height-1-m,kernel(0,m-1))); }      
      hav2TripletList.push_back(Triplet<double>(width*height-1,width*height-1-m-width-1,kernel(1,m-1)));        
     
    }
    // top/bottom edges(except corner )
    for (int j = 1; j < width-1; j++)
    {
      for (int m=0;m<kernel_size;m++ )
            {
               hav2TripletList.push_back(Triplet<double>(j,j+width-1+m,kernel(1,m)));
               if(kernel(2,m)!=0){ 
               hav2TripletList.push_back(Triplet<double>(j,j+2*width-2+m,kernel(2,m)));}     
            }
    }
    for (int j = 1; j < width-1; j++)
    {
      for (int m=0;m<kernel_size;m++ )
            {
              if(kernel(0,m)!=0)
              {
               hav2TripletList.push_back(Triplet<double>((height-2)*width+j,(height-2)*width+m,kernel(0,m)));}              
               hav2TripletList.push_back(Triplet<double>((height-1)*width+j,(height-1)*width-width-1+m,kernel(1,m))); 
              }
    }
    //other elements
    for (int i = 1; i < height-1; i++)
    {
        for (int j = 0; j < width; j++)
        {
           
            int k = i *width + j; // row number of A1
            if (j==0)
           {
            for (int m=1;m<kernel_size;++m )
              {
               
                if(kernel(0,m)!=0){
                    hav2TripletList.push_back(Triplet<double>(k,k+m,kernel(0,m)));
                }                   
                hav2TripletList.push_back(Triplet<double>(k,k-width-1+m,kernel(1,m)));               
                if(kernel(2,m)!=0){
                hav2TripletList.push_back(Triplet<double>(k,k+width-1+m,kernel(2,m))); 
                }        
              
              }
           }
           else if (j==width-1)
           {
            for (int m=0;m<kernel_size-1;++m )
              {
                if(kernel(0,m)!=0){
                    hav2TripletList.push_back(Triplet<double>(k,k+m,kernel(0,m)));
                }                   
                hav2TripletList.push_back(Triplet<double>(k,k-width-1+m,kernel(1,m)));                 
                if(kernel(2,m)!=0){
                hav2TripletList.push_back(Triplet<double>(k,k+width-1+m,kernel(2,m))); 
                }              
              }
           }
           else
           {
                for (int m=0;m<kernel_size;++m )
            {
              if(kernel(0,m)!=0)
              {
                  hav2TripletList.push_back(Triplet<double>(k,k+m,kernel(0,m)));
              }             
                  hav2TripletList.push_back(Triplet<double>(k,k-width-1+m,kernel(1,m)));                
              if(kernel(2,m)!=0)
              {
                  hav2TripletList.push_back(Triplet<double>(k,k+width-1+m,kernel(2,m)));     
              }     
               
            }
           }
            
           
        }
    }
    A.setFromTriplets(hav2TripletList.begin(), hav2TripletList.end());
    return A;
}



// Define the function that convert a vector to a Matrix<unsigned char> type and output it to image.png
void outputImage(const VectorXd &vectorData, int height, int width, const std::string &path)
{
    Matrix<double, Dynamic, Dynamic, RowMajor> output_image_matrix(height, width);
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            output_image_matrix(i, j) = vectorData(i * width + j);
        }
    }

    // Convert the modified image to grayscale and export it using stbi_write_png
    Matrix<unsigned char, Dynamic, Dynamic, RowMajor> new_image_output = output_image_matrix.unaryExpr(
        [](double pixel)
        { 
            return static_cast<unsigned char>(std::max(0.0, std::min(255.0, pixel * 255))); // ensure range [0,255]
        }
    );
    if (stbi_write_png(path.c_str(), width, height, 1, new_image_output.data(), width) == 0) {
        std::cerr << "Error: Could not save modified image" << std::endl;
    }
    std::cout << "New image saved to " << path << "\n" << std::endl;
}

// Export the vector
void exportVector(VectorXd data, const std::string &path)
{
    FILE *out = fopen(path.c_str(), "w");
    fprintf(out, "%%%%MatrixMarket vector coordinate real general\n");
    //fprintf(out, "%int\n", data.size());
    for (int i = 0; i < data.size(); i++)
    {
        fprintf(out, "%d %f\n", i, data(i));
    }
    fclose(out);
}

// Export a sparse matrix by saveMarket
void exportSparsematrix(SparseMatrix<double, RowMajor> data, const std::string &path)
{
    if (saveMarket(data, path))
    {
        std::cout << "Sparse matrix A saved to " << path << std::endl;
    }
    else
    {
        std::cerr << "Error: Could not save sparse matrix A to " << path << std::endl;
    }
}

// Check if a matrix is symmetric by tolerance 1e-10
bool isSymmetric(SparseMatrix<double, RowMajor> &matrix, const std::string &matrixName)
{
    double tolerance = 1e-10;
    double norm_diff = (matrix - SparseMatrix<double, RowMajor>(matrix.transpose())).norm();
    std::cout << matrixName << "\trow: " << matrix.rows() << "\tcolumns: " << matrix.cols() << std::endl;
    std::cout << "Check if " << matrixName << " is symmectric by norm value of its difference with transpose: "
              << norm_diff << std::endl;
    return norm_diff < tolerance;
}

// Function to check if a matrix is positive definite by cholesky
bool isPositiveDefinite(const SparseMatrix<double, RowMajor> &matrix)
{
    Eigen::SimplicialLLT<SparseMatrix<double, RowMajor>> cholesky(matrix);
    return cholesky.info() == Success;
}


/*-------------------------------------------------Main()-------------------------------------------------------*/
int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <image_path>" << std::endl;
        return 1;
    }

    const char *input_image_path = argv[1];

    /*****************************Load the image using stb_image****************************/
    int width, height, channels;
    // for greyscale images force to load only one channel
    unsigned char *image_data = stbi_load(input_image_path, &width, &height, &channels, 1);
    if (!image_data)
    {
        std::cerr << "Error: Could not load image " << input_image_path << std::endl;
        return 1;
    }

    std::cout << "Image " << argv[1] << " loaded: " << height << "x" << width << " pixels with " << channels << " channels" << std::endl;

    /****************Convert the image_data to MatrixXd form, each element value [0,1]*******************/
    // Attention! if use MatrixXd here its columnmajor by default, we prefere rowmajor
    MatrixXd noise;
    noise=MatrixXd::Random(height,width);
    noise=50*noise;
    Matrix<double, Dynamic, Dynamic, RowMajor> original_image_matrix(height, width);
    Matrix<double, Dynamic, Dynamic, RowMajor>  noised_image_matrix(height, width);
    
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            int index = (i * width + j) * channels;
            original_image_matrix(i, j) = static_cast<double>(image_data[index]) / 255;
            noised_image_matrix(i,j)=static_cast<double>(image_data[index]+noise(i,j)) / 255;
            if(noised_image_matrix(i,j)>=1) {noised_image_matrix(i,j)=1;}
            if(noised_image_matrix(i,j)<=0){noised_image_matrix(i,j)=0;}
        }
    }

    // Report the size of the matrix
    std::cout << "The original image " << argv[1] << " in matrix form has dimension: " << original_image_matrix.rows() << " rows x " << original_image_matrix.cols()
              << " cols = " << original_image_matrix.size() << "\n" << std::endl;

    /***********************************Introduce noise and export****************************************/


   

    // Convert noise image to grayscale and export it by using stbi_write_png
    Matrix<unsigned char, Dynamic, Dynamic, RowMajor> noised_image_output = noised_image_matrix.unaryExpr(
        [](double pixel)
        { return static_cast<unsigned char>(pixel * 255); });
    const std::string noised_image_path = "qiao_NoisedImage.png";
    if (stbi_write_png(noised_image_path.c_str(), width, height, 1,
                        noised_image_output.data(), width) == 0) {
        std::cerr << "Error: Could not save noised image" << std::endl;
        return 1;
    }
    std::cout << "New image saved to " << noised_image_path << "\n" << std::endl;

    /********************By map creat a vector reference to memeory without copying data***********************/
    // It is columnmajor by default, however, we've already declared our data rowmajor so here rowmajor as well.
    Map<VectorXd> v(original_image_matrix.data(), original_image_matrix.size());
    Map<VectorXd> w(noised_image_matrix.data(), noised_image_matrix.size());

    // Verify the size of the vectors
    std::cout << "Original image vector v's size: " << v.size() << std::endl;
    std::cout << "Noisy image vector w's size: " << w.size() << std::endl;
    std::cout << "Euclidean norm of v is: " << v.norm() << std::endl;

    /****************************Create different kernels and perform convolution*****************************/
    const int kernel_size = 3;

    // Create the kernel H_{av2}
    Matrix<double, kernel_size, kernel_size, RowMajor> hav2;
    hav2 << 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0,
            1.0, 1.0, 1.0;
    hav2=(1/9.0)*hav2;
    SparseMatrix<double, RowMajor> A1 = convolution3(hav2, height, width,3);
    // Check the nonzero numbers
    std::cout << "A1 nonzero number is: " << A1.nonZeros() << std::endl;
    // Smooth the noise image by using this filterImage function and passing data and path into it
    const std::string smooth_image_path = "qiao_SmoothedImage.png";
    VectorXd smoothed_image_vector = A1 * w;
    outputImage(smoothed_image_vector, height, width, smooth_image_path);

    // Create the kernel H_{sh2}}
    Matrix<double, kernel_size, kernel_size, RowMajor> hsh2;
    hsh2 << 0.0, -3.0, 0.0,
           -1.0, 9.0, -3.0,
            0.0, -1.0, 0.0;
    // Perform convolution function
    SparseMatrix<double, RowMajor> A2 = convolution3(hsh2, height, width,3);
    // Check the nonzero numbers
    std::cout << "A2 nonzero number is: " << A2.nonZeros() << std::endl;
    // Verify if the matrix A2 is symmetric
    isSymmetric(A2, "A2") ? std::cout << "The matrix A2 is symmetric!" << std::endl
                          : std::cout << "The matrix A2 is not symmetric!" << std::endl;

    // Sharpen the original image by matrix production
    const std::string sharpen_image_path = "qiao_SharpenedImage.png";
    VectorXd sharpened_image_vector = A2 * v;
    outputImage(sharpened_image_vector, height, width, sharpen_image_path);

    // Create the kernel H_{sh2}}
    Matrix<double, kernel_size, kernel_size, RowMajor> hlap;
    hlap << 0.0, -1.0, 0.0,
           -1.0, 4.0, -1.0,
            0.0, -1.0, 0.0;
    // Perform convolution function
    SparseMatrix<double, RowMajor> A3 = convolution3(hlap, height, width,3);
    // Verify if the matrix A3 is symmetric
    isSymmetric(A3, "A3") ? std::cout << "The matrix A3 is symmetric!" << std::endl
                          : std::cout << "The matrix A3 is not symmetric!" << std::endl;
    
    // Check if A3 is positive definite as well
    isPositiveDefinite(A3) ? std::cout << "The matrix A3 is positive definite!" << std::endl
                           : std::cout << "The matrix A3 is not positive definite!" << std::endl;

    // Edge detection of the original image
    const std::string edgeDetection_image_path = "qiao_EdgeDetectionImage.png";
    VectorXd edgeDetected_sharpened_image_vector = A3 * v;
    outputImage(edgeDetected_sharpened_image_vector, height, width, edgeDetection_image_path);

    
     VectorXd mat(height*width);
     VectorXd mat2(height*width);
     loadMarketVector(mat,"pre_x.mtx");
     
     
     /**********************************Solve equation of A2*x = w****************************************/
    
    // Read x.mtx file data and output it as image
    const std::string x_image_path = "output_VectorX.png";
    outputImage(mat, height, width, x_image_path);
    return 0;
}