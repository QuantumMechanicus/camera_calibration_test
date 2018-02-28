#include <Eigen/Dense>
#include <tbb/tbb.h>
#include <opencv2/core/core.hpp>
#include <opencv2/core/eigen.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <boost/program_options.hpp>
#include "Core.h"


int main(int argc, char *argv[]) {
    std::string input_image_name;
    std::vector<double> v_lambdas;
    std::string out_path;
    namespace po = boost::program_options;
    try {
        // clang-format off
        po::options_description desc("Global reconstruction options");
        desc.add_options()
                ("help", "Print help message")
                ("i", po::value<std::string>(&input_image_name)->required(),
                 "Input image name")
                ("o", po::value<std::string>(&out_path)->default_value("./Undistorted.JPG"),
                 "Output image name")
                ("c", po::value<std::vector<double> >(&v_lambdas)->multitoken()->required(), "Distortion coefficients");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc,
                                         po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);

        // clang-format on
        if (vm.count("help")) {
            std::cout << desc << std::endl;
            return -1;
        }
        boost::program_options::notify(vm);

    } catch (std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -2;
    }

    cv::Mat in, out, mx, my;
    in = cv::imread(input_image_name);
    int rows = in.rows;
    int cols = in.cols;


    double r_img = std::sqrt(cols * cols / 4.0 + rows * rows / 4.0);

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> map_x(rows, cols), map_y(rows, cols);

    double d = std::max(rows, cols) / 2.0;
    double dr = d / r_img;

    size_t first_non_zero = v_lambdas.size() - 1;
    while (v_lambdas[first_non_zero] == 0)
        --first_non_zero;
    v_lambdas.resize(first_non_zero + 1);
    Eigen::VectorXd lambdas = Eigen::Map<Eigen::VectorXd>(v_lambdas.data(), v_lambdas.size());


    auto alpha = utils::distortion_problem::undistortionDenominator<double>(dr, lambdas.cast<double>());


    tbb::parallel_for(tbb::blocked_range<int>(0, cols), [&](auto range) {

        for (int i = range.begin(); i != range.end(); ++i) {
            for (int j = 0; j < rows; ++j) {
                double ii = ((i - cols / 2.0) / r_img) / alpha;
                double jj = ((j - rows / 2.0) / r_img) / alpha;
                double r = std::sqrt(ii * ii + jj * jj + 1e-7);
                double rd = utils::distortion_problem::findDistortedRadius(lambdas, r);
                auto dd = utils::distortion_problem::undistortionDenominator<double>(rd, lambdas.cast<double>());
                map_x(j, i) = r_img * ii * dd + cols / 2.0;
                map_y(j, i) = r_img * jj * dd + rows / 2.0;
            }
        }

    });

    cv::eigen2cv(map_x.cast<float>().eval(), mx);
    cv::eigen2cv(map_y.cast<float>().eval(), my);
    cv::remap(in, out, mx, my, cv::INTER_LINEAR, cv::BORDER_CONSTANT, cv::Scalar(0, 0, 0));
    cv::imwrite(out_path, out);
    std::cout << "Undistortion've done\n";

    return 0;
}
