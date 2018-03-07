//
// Created by danielbord on 2/14/18.
//
#include <boost/program_options.hpp>
#include "Global_Non_Linear_Estimator.h"

int main(int argc, char *argv[]) {

    std::vector<std::string> f_infos;
    int number_of_distortion_coefficients;
    int non_linear_iter;
    double percent_of_inliers;

    namespace po = boost::program_options;
    try {

        po::options_description desc("Automatic solver options");
        desc.add_options()
                ("help", "Print help message")
                ("d", po::value<std::vector<std::string>>(&f_infos)->multitoken()->required(),
                 "File with data describes path to keypoints, cameras info etc (see examples)")
                ("i", po::value<int>(&non_linear_iter)->default_value(1), "Number of non linear iterations")
                ("q", po::value<double>(&percent_of_inliers)->default_value(0.1), "quantile to minimize")
                ("c", po::value<int>(&number_of_distortion_coefficients)->default_value(1),
                 "number of distortion coefficients");

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

    std::shared_ptr<scene::DynamicDivisionModelStereoPair::VertexMap_t> cameras
            = std::make_shared<scene::DynamicDivisionModelStereoPair::VertexMap_t>(
                    scene::DynamicDivisionModelStereoPair::VertexMap_t());

    std::vector<scene::DynamicDivisionModelStereoPair> stereo_pairs(number_of_pairs);

    double mean_distortion_coefficient = 0;
    double mean_ppx = 0;
    double mean_ppy = 0;
    double mean_f = 0;
    double w = 0;
    double h = 0;

    std::vector<scene_serialization::DynamicSceneArchiver> archivers(
            number_of_pairs);
    for (size_t k = 0; k < number_of_pairs; ++k) {
        archivers[k].parse(f_infos[k]);
        scene::DynamicDivisionModelStereoPair stereo_pair;
        stereo_pair.loadScene(archivers[k], cameras);
        assert(stereo_pair.getLeftIntrinsicsPointer()->getNumberOfCoefficients() > 0 &&
               stereo_pair.getRightIntrinsicsPointer()->getNumberOfCoefficients() > 0 && "Incorrect intrinsics data");
        stereo_pairs[k] = stereo_pair;
        mean_distortion_coefficient += (stereo_pair.getLeftIntrinsics().getDistortionCoefficients()(0) +
                                        stereo_pair.getRightIntrinsics().getDistortionCoefficients()(0)) / 2.0;
        mean_ppx += (stereo_pair.getLeftIntrinsics().getPrincipalPointX() +
                     stereo_pair.getRightIntrinsics().getPrincipalPointX()) / 2.0;
        mean_ppy += (stereo_pair.getLeftIntrinsics().getPrincipalPointY() +
                     stereo_pair.getRightIntrinsics().getPrincipalPointY()) / 2.0;
        mean_f += (stereo_pair.getLeftIntrinsics().getFocalLength() +
                   stereo_pair.getRightIntrinsics().getFocalLength()) / 2.0;

        w += (stereo_pairs[k].getLeftIntrinsicsPointer()->getWidth() +
              stereo_pairs[k].getRightIntrinsicsPointer()->getWidth()) / 2.0;

        h += (stereo_pairs[k].getLeftIntrinsicsPointer()->getHeight() +
              stereo_pairs[k].getRightIntrinsicsPointer()->getHeight()) / 2.0;

        stereo_pairs[k].normalizeLeftKeypoints();
        stereo_pairs[k].normalizeRightKeypoints();

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
    std::shared_ptr<intrinsics::DynamicDivisionModel> common_intrinsics_parameters = std::make_shared<intrinsics::DynamicDivisionModel>(
            distortion_coefficients,
            w, h, mean_f, mean_ppx, mean_ppy);

    scene::StdVector<scene::FundamentalMatrix> fundamental_matrices(number_of_pairs);
    scene::StdVector<Eigen::Vector3d> translations(number_of_pairs);
    scene::StdVector<Eigen::Matrix3d> rotations(number_of_pairs);
    scene::StdVector<scene::ImagePoints> left_pictures_keypoints(number_of_pairs), right_pictures_keypoints(
            number_of_pairs);

    for (size_t k = 0; k < number_of_pairs; ++k) {

        left_pictures_keypoints[k] = stereo_pairs[k].getLeftKeypoints();
        right_pictures_keypoints[k] = stereo_pairs[k].getRightKeypoints();
        fundamental_matrices[k] = stereo_pairs[k].getFundamentalMatrix();
        stereo_pairs[k].recoverRelativeMotion();
        translations[k] = stereo_pairs[k].getRelativeTranslation();
        rotations[k] = stereo_pairs[k].getRelativeRotation().matrix();
        stereo_pairs[k].estimateLeftCamera(common_intrinsics_parameters);
        stereo_pairs[k].estimateRightCamera(common_intrinsics_parameters);
    }
    double image_radius = std::sqrt((w / 2) * (w / 2) + (h / 2) * (h / 2));
    non_linear_optimization::GlobalNonLinearEstimatorOptions options(percent_of_inliers, 100, image_radius);
    non_linear_optimization::GlobalNonLinearEstimator estimator(left_pictures_keypoints, right_pictures_keypoints,
                                                                common_intrinsics_parameters->getDistortionCoefficients(),
                                                                rotations, translations,
                                                                common_intrinsics_parameters->getFocalLength(),
                                                                common_intrinsics_parameters->getPrincipalPointX(),
                                                                common_intrinsics_parameters->getPrincipalPointY(),
                                                                options);


    common_intrinsics_parameters->estimateParameter(estimator);
    scene::DynamicDivisionModelScene scene(cameras, stereo_pairs);
    scene.estimateStereoPairs(estimator);
    scene.saveScene(archivers);

    return 0;
}