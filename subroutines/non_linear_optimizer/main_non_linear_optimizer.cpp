//
// Created by danielbord on 2/14/18.
//
#include <boost/program_options.hpp>
#include "Non_Linear_Estimator.h"

int main(int argc, char *argv[]) {

    std::vector<std::string> f_infos;
    int non_linear_iter;
    double percent_of_inliers;
    int number_of_distortion_coefficients;
    namespace po = boost::program_options;
    try {

        po::options_description desc("Automatic solver options");
        desc.add_options()
                ("help", "Print help message")
                ("d", po::value<std::vector<std::string>>(&f_infos)->multitoken()->required(),
                 "File with data describes path to keypoints, cameras info etc (see examples)")
                ("i", po::value<int>(&non_linear_iter)->default_value(1), "Number of non linear iterations")
                ("q", po::value<double>(&percent_of_inliers)->default_value(0.1), "quantile to minimize"),
                ("c", po::value<int>(&number_of_distortion_coefficients)->default_value(1), "number of distortion coefficients");

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
    auto number_of_pairs = f_infos.size();


    std::shared_ptr<std::map<std::string, scene::Camera<intrinsics::DivisionModelIntrinsic<Eigen::Dynamic>>>> cameras
            = std::make_shared<std::map<std::string, scene::Camera<intrinsics::DivisionModelIntrinsic<Eigen::Dynamic>>>>(
                    std::map<std::string, scene::Camera<intrinsics::DivisionModelIntrinsic<Eigen::Dynamic>>>());

    std::vector<scene::TwoView<intrinsics::DivisionModelIntrinsic<Eigen::Dynamic>>> stereo_pairs(number_of_pairs);

    double mean_distortion_coefficient = 0;
    double mean_ppx = 0;
    double mean_ppy = 0;
    double mean_f = 0;
    double w = 0;
    double h = 0;
    for (size_t k = 0; k < number_of_pairs; ++k) {
        scene_serialization::SimpleSceneArchiver<scene_serialization::SimpleDivisionModelArchiver<std::string, -1>> archiver;
        archiver.parse(f_infos[k]);
        scene::TwoView<intrinsics::DivisionModelIntrinsic<Eigen::Dynamic>> stereo_pair;
        stereo_pair.loadScene(archiver, cameras);
        stereo_pairs[k] = stereo_pair;
        mean_distortion_coefficient += stereo_pair.getLeftIntrinsics().getDistortionCoefficients()(0);
        mean_ppx += stereo_pair.getLeftIntrinsics().getPrincipalPointX();
        mean_ppy += stereo_pair.getLeftIntrinsics().getPrincipalPointY();
        mean_f += stereo_pair.getLeftIntrinsics().getFocalLength();
        stereo_pairs[k].normalizeLeftKeypoints();
        stereo_pairs[k].normalizeRightKeypoints();
        w += (stereo_pairs[k].getLeftIntrinsicsPointer()->getWidth() +
              stereo_pairs[k].getLeftIntrinsicsPointer()->getWidth()) / 2.0;

        h += (stereo_pairs[k].getLeftIntrinsicsPointer()->getHeight() +
              stereo_pairs[k].getLeftIntrinsicsPointer()->getHeight()) / 2.0;

    }
    w /= number_of_pairs;
    h /= number_of_pairs;

    mean_distortion_coefficient /= number_of_pairs;
    mean_ppx /= number_of_pairs;
    mean_ppy /= number_of_pairs;
    mean_f /= number_of_pairs;
    Eigen::RowVectorXd distortion_coefficients(number_of_distortion_coefficients);
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

        stereo_pairs[k].estimateLeftIntrinsics(common_intrinsics_parameters);
        stereo_pairs[k].estimateRightIntrinsics(common_intrinsics_parameters);
    }
    double image_radius = std::sqrt((w / 2) * (w / 2) + (h / 2) * (h / 2));
    non_linear_optimization::NonLinearEstimatorOptions options(non_linear_iter, percent_of_inliers, 100, image_radius);
    non_linear_optimization::NonLinearEstimator estimator(left_pictures_keypoints, right_pictures_keypoints,
                                                          fundamental_matrices,
                                                          common_intrinsics_parameters->getDistortionCoefficients(),
                                                          options);
    common_intrinsics_parameters->estimateParameter(estimator);
    scene::Scene<scene::Camera<intrinsics::DivisionModelIntrinsic<-1>>> scene(cameras, stereo_pairs);
    scene.estimateStereoPairs(estimator);
    /*
    //This is for test (later will be added for gtest)
    std::cout << "T\n";
    std::cout << common_intrinsics_parameters.use_count() << std::endl;
    for (size_t k = 0; k < number_of_pairs; ++k) {

        stereo_pairs[k] = scene::TwoView<intrinsics::DivisionModelIntrinsic<Eigen::Dynamic>>(
                scene::Camera<intrinsics::DivisionModelIntrinsic<Eigen::Dynamic>>(common_intrinsics_parameters,
                                                                                           stereo_pairs[k].getLeftRotation(),
                                                                                           stereo_pairs[k].getLeftTranslation()),
                scene::Camera<intrinsics::DivisionModelIntrinsic<Eigen::Dynamic>>(common_intrinsics_parameters,
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