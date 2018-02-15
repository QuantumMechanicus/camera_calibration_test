//
// Created by danielbord on 2/14/18.
//
#include <boost/program_options.hpp>
#include "Non_Linear_Estimator.h"

int main(int argc, char *argv[]) {

    unsigned int w, h;
    std::vector<std::string> input1, input2, f_estimated_fundamental_matrices, f_estimated_intrinsics;
    int non_linear_iter, number_of_distortion_coefficients;
    double percent_of_inliers;
    std::string f_optimizer_results;
    namespace po = boost::program_options;
    try {

        po::options_description desc("Optimization input options");
        desc.add_options()
                ("help", "Print help message")
                ("ff",
                 po::value<std::vector<std::string> >(&f_estimated_fundamental_matrices)->multitoken()->required(),
                 "Path to estimated fundamental matrices")
                ("lf",
                 po::value<std::vector<std::string> >(&f_estimated_intrinsics)->multitoken()->required(),
                 "Path to estimated intrinsics")
                ("lp", po::value<std::vector<std::string> >(&input1)->multitoken()->required(),
                 "Path to files with first (left) camera keypoints")
                ("rp", po::value<std::vector<std::string> >(&input2)->multitoken()->required(),
                 "Path to files with second (right) camera keypoints")
                ("rf", po::value<std::string>(&f_optimizer_results)->default_value("./optimizer_results/"),
                 "Path to results directory")
                ("i", po::value<int>(&non_linear_iter)->default_value(1), "Number of iterations")
                ("w", po::value<unsigned int>(&w)->required(), "Width")
                ("h", po::value<unsigned int>(&h)->required(), "Height")
                ("q", po::value<double>(&percent_of_inliers)->default_value(0.1), "quantile to minimize")
                ("nl", po::value<int>(&number_of_distortion_coefficients)->default_value(1),
                 "Number of parameters in denominator of model");

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

    /*assert(input1.size() == input2.size() && "number of left and right cameras should be equal");
    assert(f_estimated_fundamental_matrices.size() == input1.size() &&
           "number of cameras and fundamental matrices should be equal");*/
    long number_of_pairs = f_estimated_fundamental_matrices.size();


    std::vector<scene::TwoView<intrinsics::DivisionModelIntrinsic<Eigen::Dynamic>>> stereo_pairs(number_of_pairs);
    double mean_distortion_coefficient = 0;
    double mean_ppx = 0;
    double mean_ppy = 0;
    double mean_f = 0;
    for (size_t k = 0; k < number_of_pairs; ++k) {
        scene_serialization::SimpleSceneArchiver<scene_serialization::SimpleDivisionModelArchiver<Eigen::Dynamic>> archiver(
                f_estimated_fundamental_matrices[k],
                f_estimated_intrinsics[k], f_estimated_intrinsics[k], input1[k], input2[k]);
        scene::TwoView<intrinsics::DivisionModelIntrinsic<Eigen::Dynamic>> stereo_pair;
        stereo_pair.loadScene(archiver);
        stereo_pairs[k] = stereo_pair;
        mean_distortion_coefficient += stereo_pair.getLeftIntrinsics().getDistortionCoefficients()(0);
        mean_ppx += stereo_pair.getLeftIntrinsics().getPrincipalPointX();
        mean_ppy += stereo_pair.getLeftIntrinsics().getPrincipalPointY();
        mean_f += stereo_pair.getLeftIntrinsics().getFocalLength();
        stereo_pairs[k].normalizeLeftKeypoints();
        stereo_pairs[k].normalizeRightKeypoints();

    }


    mean_distortion_coefficient /= number_of_pairs;
    mean_ppx /= number_of_pairs;
    mean_ppy /= number_of_pairs;
    mean_f /= number_of_pairs;
    Eigen::VectorXd distortion_coefficients(number_of_distortion_coefficients);
    distortion_coefficients.setZero();
    distortion_coefficients[0] = mean_distortion_coefficient;
    std::shared_ptr<intrinsics::DivisionModelIntrinsic<Eigen::Dynamic>> common_intrinsics_parameters = std::make_shared<intrinsics::DivisionModelIntrinsic<Eigen::Dynamic>>(
            distortion_coefficients,
            w, h, mean_f, mean_ppx, mean_ppy);

    scene::StdVector<scene::FundamentalMatrix> fundamental_matrices;
    scene::StdVector<scene::ImagePoints> left_pictures_keypoints, right_pictures_keypoints;

    for (size_t k = 0; k < number_of_pairs; ++k) {
        left_pictures_keypoints.push_back(stereo_pairs[k].getLeftKeypoints());
        right_pictures_keypoints.push_back(stereo_pairs[k].getRightKeypoints());
        fundamental_matrices.push_back(stereo_pairs[k].getFundamentalMatrix());

        stereo_pairs[k] = scene::TwoView<intrinsics::DivisionModelIntrinsic<Eigen::Dynamic>>(
                internal_scene::Camera<intrinsics::DivisionModelIntrinsic<Eigen::Dynamic>>(common_intrinsics_parameters,
                                                                                           stereo_pairs[k].getLeftRotation(),
                                                                                           stereo_pairs[k].getLeftTranslation()),
                internal_scene::Camera<intrinsics::DivisionModelIntrinsic<Eigen::Dynamic>>(common_intrinsics_parameters,
                                                                                           stereo_pairs[k].getRightRotation(),
                                                                                           stereo_pairs[k].getRightTranslation()),
                stereo_pairs[k].getLeftKeypoints(), stereo_pairs[k].getRightKeypoints(),
                stereo_pairs[k].getFundamentalMatrix());
    }
    double image_radius = std::sqrt((w / 2) * (w / 2) + (h / 2) * (h / 2));
    non_linear_optimization::NonLinearEstimatorOptions options(non_linear_iter, percent_of_inliers, 100, image_radius);
    non_linear_optimization::NonLinearEstimator estimator(left_pictures_keypoints, right_pictures_keypoints,
                                                          fundamental_matrices, *common_intrinsics_parameters, options);
    estimator.estimate();
    /*
    //This is for test (later will be added for gtest)
    std::cout << "T\n";
    std::cout << common_intrinsics_parameters.use_count() << std::endl;
    for (size_t k = 0; k < number_of_pairs; ++k) {

        stereo_pairs[k] = scene::TwoView<intrinsics::DivisionModelIntrinsic<Eigen::Dynamic>>(
                internal_scene::Camera<intrinsics::DivisionModelIntrinsic<Eigen::Dynamic>>(common_intrinsics_parameters,
                                                                                           stereo_pairs[k].getLeftRotation(),
                                                                                           stereo_pairs[k].getLeftTranslation()),
                internal_scene::Camera<intrinsics::DivisionModelIntrinsic<Eigen::Dynamic>>(common_intrinsics_parameters,
                                                                                           stereo_pairs[k].getRightRotation(),
                                                                                           stereo_pairs[k].getRightTranslation()),
                stereo_pairs[k].getLeftKeypoints(), stereo_pairs[k].getRightKeypoints());
        std::cout << stereo_pairs[k].getLeftIntrinsicsPointer().use_count() << std::endl;
        std::cout << stereo_pairs[k].getRightIntrinsicsPointer().use_count() << std::endl;

    }
    common_intrinsics_parameters.reset();
    std::cout << stereo_pairs[0].getLeftIntrinsicsPointer().use_count() << std::endl;
    std::cout << stereo_pairs[0].getRightIntrinsicsPointer().use_count() << std::endl;

    /*scene_serialization::SimpleSceneArchiver<scene_serialization::SimpleDivisionModelArchiver<Eigen::Dynamic>>
            archiver("/home/danielbord/CLionProjects/camera_calibration/pipeline/bin_debug/estimated_f",
                     "/home/danielbord/CLionProjects/camera_calibration/pipeline/bin_debug/estimated_lambda","/home/danielbord/CLionProjects/camera_calibration/pipeline/bin_debug/estimated_lambda");

    scene::TwoView<intrinsics::DivisionModelIntrinsic<Eigen::Dynamic>> stereo_pair2;
    stereo_pair2.loadScene(archiver);
    std::cout << stereo_pair2.getFundamentalMatrix() << std::endl;
    std::cout << stereo_pair2.getLeftIntrinsics().getHeight() << std::endl;

    std::cout << stereo_pair2.getLeftIntrinsics().getDistortionCoefficients() << std::endl;*/
    return 0;
}