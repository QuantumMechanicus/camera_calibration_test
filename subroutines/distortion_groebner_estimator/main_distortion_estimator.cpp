#include <chrono>
#include <ctime>
#include <fstream>
#include <iostream>
#include <boost/program_options.hpp>
#include "Groebner_Estimator.h"


int main(int argc, char *argv[]) {
    double w, h, r;
    std::string input1, input2, f_estimated_lambda, f_estimated_fundamental_matrix, f_inliers_list;
    int iter;
    double lambda_upper_bound;
    double lambda_lower_bound;
    double percent_of_inliers;
    namespace po = boost::program_options;
    try {

        po::options_description desc("Automatic solver options");
        desc.add_options()
                ("help", "Print help message")
                ("lp", po::value<std::string>(&input1)->required(),
                 "Path to file with first (left) camera keypoints")
                ("rp", po::value<std::string>(&input2)->required(),
                 "Path to file with second (right) camera keypoints")
                ("up_threshold", po::value<double>(&lambda_upper_bound)->default_value(0.25), "Lambda upper threshold")
                ("low_threshold", po::value<double>(&lambda_lower_bound)->default_value(-1), "Lambda lower threshold")
                ("ff", po::value<std::string>(&f_estimated_fundamental_matrix)->default_value(
                        "./estimated_f"), "Output file for fundamental matrix estimation")
                ("lf", po::value<std::string>(&f_estimated_lambda)->default_value(
                        "./estimated_lambda"), "Output file for lambda estimation")
                ("i", po::value<int>(&iter)->default_value(10000), "Number of iterations")
                ("w", po::value<double>(&w)->default_value(7360), "Width")
                ("h", po::value<double>(&h)->default_value(4912), "Height")
                ("q", po::value<double>(&percent_of_inliers)->default_value(0.1), "quantile to minimize")
                ("if", po::value<std::string>(&f_inliers_list)->default_value("./automatic_solver_results/inliers"),
                 "Output file for inliers");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);

        if (vm.count("help")) {
            std::cout << desc << std::endl;
            return -1;
        }
        boost::program_options::notify(vm);

    } catch (std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -2;
    }
    r = std::sqrt(w * w + h * h) / 2;

    Eigen::Matrix<double, 2, Eigen::Dynamic> u1d, u2d;
    utils::loadMatrix<double, 2, Eigen::Dynamic>(input1, u1d, true);
    utils::loadMatrix<double, 2, Eigen::Dynamic>(input2, u2d, true);

    Eigen::Matrix3d fundamental_matrix;
    std::shared_ptr<intrinsics::StandardDivisionModelIntrinsic> intrinsics = std::make_shared<intrinsics::StandardDivisionModelIntrinsic>(
            w, h);

    internal_scene::Camera<intrinsics::StandardDivisionModelIntrinsic> left_camera(intrinsics);
    internal_scene::Camera<intrinsics::StandardDivisionModelIntrinsic> right_camera(intrinsics);

    scene::TwoView<intrinsics::DivisionModelIntrinsic<1>> stereo_pair(left_camera, right_camera, u1d, u2d);
    stereo_pair.normalizeLeftKeypoints();
    stereo_pair.normalizeRightKeypoints();

    auto start = std::chrono::high_resolution_clock::now();

    estimators::GroebnerEstimatorOptions opt(iter, percent_of_inliers, lambda_lower_bound, lambda_upper_bound, r);
    estimators::GroebnerDivisionModelEstimator groebner_estimator(stereo_pair.getLeftKeypoints(),
                                                                  stereo_pair.getRightKeypoints(), opt);
    stereo_pair.estimateLeftIntrinsics(groebner_estimator);
    stereo_pair.estimateFundamentalMatrix(groebner_estimator);
    scene_serialization::SimpleSceneArchiver<scene_serialization::SimpleDivisionModelArchiver<1>>
            archiver(f_estimated_fundamental_matrix, f_estimated_lambda, f_estimated_lambda);
    stereo_pair.saveScene(archiver);


    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start).count();

    std::cout << std::endl << "Estimation done in " << std::setprecision(5) << duration / 1000.0 << " seconds"
              << std::endl;

    //TODO add inliers list
    return 0;
}