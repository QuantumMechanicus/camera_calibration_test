#include <chrono>
#include <ctime>
#include <fstream>
#include <iostream>
#include <boost/program_options.hpp>
#include "Groebner_Estimator.h"


int main(int argc, char *argv[]) {
    std::string f_inliers_list, f_info;
    int iter;
    double lambda_upper_bound;
    double lambda_lower_bound;
    double percent_of_inliers;
    namespace po = boost::program_options;
    try {

        po::options_description desc("Automatic solver options");
        desc.add_options()
                ("help", "Print help message")
                ("up_threshold", po::value<double>(&lambda_upper_bound)->default_value(0.25), "Lambda upper threshold")
                ("low_threshold", po::value<double>(&lambda_lower_bound)->default_value(-1), "Lambda lower threshold")
                ("d", po::value<std::string>(&f_info)->required(),
                 "File with data describes path to keypoints, cameras info etc (see examples)")
                ("i", po::value<int>(&iter)->default_value(10000), "Number of iterations")
                ("q", po::value<double>(&percent_of_inliers)->default_value(0.1), "quantile to minimize")
                ("if",
                 po::value<std::string>(&f_inliers_list)->default_value("./automatic_solver_results/inliers"),
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


    Eigen::Matrix3d fundamental_matrix;
    scene::TwoView<intrinsics::DivisionModelIntrinsic<>> stereo_pair;
    scene_serialization::SimpleSceneArchiver<scene_serialization::SimpleDivisionModelArchiver<>> archiver;
    archiver.parse(f_info);
    stereo_pair.loadScene(archiver);
    double w =
            (stereo_pair.getLeftIntrinsicsPointer()->getWidth() + stereo_pair.getRightIntrinsicsPointer()->getWidth()) /
            2;
    double h = (stereo_pair.getLeftIntrinsicsPointer()->getHeight() +
                stereo_pair.getRightIntrinsicsPointer()->getHeight()) / 2;
    double r = std::sqrt(w * w + h * h) / 2;
    if (w > 0 and h > 0) {
        stereo_pair.normalizeLeftKeypoints();
        stereo_pair.normalizeRightKeypoints();

        auto start = std::chrono::high_resolution_clock::now();

        estimators::GroebnerEstimatorOptions opt(iter, percent_of_inliers, lambda_lower_bound, lambda_upper_bound, r);
        estimators::GroebnerDivisionModelEstimator groebner_estimator(stereo_pair.getLeftKeypoints(),
                                                                      stereo_pair.getRightKeypoints(), opt);
        stereo_pair.estimateLeftIntrinsics(groebner_estimator);
        stereo_pair.estimateRightIntrinsics(groebner_estimator);
        stereo_pair.estimateFundamentalMatrix(groebner_estimator);


        stereo_pair.denormalizeLeftKeypoints();
        stereo_pair.denormalizeRightKeypoints();
        stereo_pair.saveScene(archiver);


        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - start).count();

        std::cout << std::endl << "Estimation done in " << std::setprecision(5) << duration / 1000.0 << " seconds"
                  << std::endl;

        //TODO add inliers list
    } else {

        std::cerr
                << "You should specify width, height and initial model parameters in files described in info (data) file (see examples)\n";
    }
    return 0;
}