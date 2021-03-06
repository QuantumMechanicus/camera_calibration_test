//
// Created by danielbord on 2/14/18.
//
#include <boost/program_options.hpp>
#include "Global_Non_Linear_Estimator.h"

#include <glog/logging.h>

int main(int argc, char *argv[]) {
    google::InitGoogleLogging(argv[0]);
    google::InstallFailureSignalHandler();
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
                ("c", po::value<int>(&number_of_distortion_coefficients)->default_value(2),
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

    auto cameras
            = std::make_shared<scene::DynamicDivisionModelStereoPair::VertexMap_t>();

    std::vector<scene::DynamicDivisionModelStereoPair> stereo_pairs(number_of_pairs);

    Eigen::VectorXd mean_distortion_coefficient;
    double mean_ppx = 0;
    double mean_ppy = 0;
    double mean_f = 0;
    double w = 0;
    double h = 0;

    std::vector<scene_serialization::DynamicSceneArchiver> archivers(
            number_of_pairs);
    for (size_t k = 0; k < number_of_pairs; ++k) {
        archivers[k].parse(f_infos[k]); //
        scene::DynamicDivisionModelStereoPair stereo_pair;
        stereo_pair.loadScene(archivers[k], cameras);
        assert(stereo_pair.getLeftIntrinsicsPointer()->getNumberOfCoefficients() > 0 &&
               stereo_pair.getRightIntrinsicsPointer()->getNumberOfCoefficients() > 0 && "Incorrect intrinsics data");
        stereo_pairs[k] = stereo_pair;
        if (k == 0) {
            mean_distortion_coefficient = Eigen::VectorXd::Zero(
                    stereo_pair.getLeftIntrinsicsPointer()->getDistortionCoefficients().rows());
        }
        mean_distortion_coefficient += (stereo_pair.getLeftIntrinsics().getDistortionCoefficients() +
                                        stereo_pair.getRightIntrinsics().getDistortionCoefficients()) / 2.0;
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


    Eigen::VectorXd distortion_coefficients(number_of_distortion_coefficients);
    distortion_coefficients.setZero();
    distortion_coefficients.head(
            std::min(static_cast<int>(mean_distortion_coefficient.rows()),
                     number_of_distortion_coefficients)) = mean_distortion_coefficient.head(
            std::min(static_cast<int>(mean_distortion_coefficient.rows()), number_of_distortion_coefficients));
    std::shared_ptr<intrinsics::DynamicDivisionModel> common_intrinsics_parameters = std::make_shared<intrinsics::DynamicDivisionModel>(
            distortion_coefficients,
            w, h, mean_f, 0, 0);

    scene::StdVector<scene::FundamentalMatrix> fundamental_matrices(number_of_pairs);
    scene::StdVector<Eigen::Vector3d> translations(number_of_pairs);
    scene::StdVector<Sophus::SO3d> rotations(number_of_pairs);
    scene::StdVector<scene::ImagePoints> left_pictures_keypoints(number_of_pairs), right_pictures_keypoints(
            number_of_pairs);
    double image_radius = std::sqrt((w / 2) * (w / 2) + (h / 2) * (h / 2));
    Eigen::Matrix3d calibration_matrix = common_intrinsics_parameters->getCalibrationMatrix();

    std::cout << calibration_matrix << std::endl;
    for (size_t k = 0; k < number_of_pairs; ++k) {
        stereo_pairs[k].estimateLeftCamera(common_intrinsics_parameters);
        stereo_pairs[k].estimateRightCamera(common_intrinsics_parameters);
        left_pictures_keypoints[k] = stereo_pairs[k].getLeftKeypoints();
        right_pictures_keypoints[k] = stereo_pairs[k].getRightKeypoints();
        fundamental_matrices[k] = stereo_pairs[k].getFundamentalMatrix();
        std::vector<size_t> inliers_ind;
        double interval = utils::distortion_problem::findInliers(stereo_pairs[k].getLeftKeypoints(),
                                                                 stereo_pairs[k].getRightKeypoints(),
                                                                 stereo_pairs[k].getLeftIntrinsicsPointer()->getDistortionCoefficients(),
                                                                 fundamental_matrices[k],
                                                                 0.1,
                                                                 inliers_ind, image_radius);
        std::cout << "Interval: " << interval << " " << image_radius << std::endl;

        //Eigen::Matrix3d hhh = stereo_pairs[k].getEssentialMatrix();

        stereo_pairs[k].recoverRelativeMotion();

        translations[k] = stereo_pairs[k].getRelativeTranslation();
        rotations[k] = stereo_pairs[k].getRelativeRotation();

        std::cout << fundamental_matrices[k] << std::endl;
        Eigen::Matrix3d fff =
                stereo_pairs[k].getRightIntrinsicsPointer()->getCalibrationMatrix().transpose().inverse() *
                utils::screw_hat(translations[k]) * rotations[k].matrix() *
                stereo_pairs[k].getLeftIntrinsicsPointer()->getCalibrationMatrix().inverse();

        std::cout << fff / fff(2, 2) << std::endl;
        std::cout << fundamental_matrices[k] << std::endl;
        Sophus::SE3d mm(stereo_pairs[k].getRelativeRotation(), stereo_pairs[k].getRelativeTranslation());
        double interval2 = utils::distortion_problem::findInliers(stereo_pairs[k].getLeftKeypoints(),
                                                                  stereo_pairs[k].getRightKeypoints(),
                                                                  stereo_pairs[k].getLeftIntrinsicsPointer()->getDistortionCoefficients(),
                                                                  mm,
                                                                  stereo_pairs[k].getLeftIntrinsicsPointer()->getCalibrationMatrix(),
                                                                  0.1,
                                                                  inliers_ind, image_radius);
        std::cout << "Interval2: " << interval2 << " " << image_radius << std::endl;
    }
    std::cout << "FCL: " << common_intrinsics_parameters->getFocalLength() << std::endl;

    non_linear_optimization::GlobalNonLinearEstimatorOptions options(percent_of_inliers, 100, image_radius, common_intrinsics_parameters->getWidth(),
    common_intrinsics_parameters->getHeight());
    non_linear_optimization::GlobalNonLinearEstimator estimator(left_pictures_keypoints, right_pictures_keypoints,
                                                                common_intrinsics_parameters->getDistortionCoefficients(),
                                                                rotations, translations,
                                                                common_intrinsics_parameters->getFocalLength(),
                                                                common_intrinsics_parameters->getPrincipalPointX(),
                                                                common_intrinsics_parameters->getPrincipalPointY(),
                                                                options);


    common_intrinsics_parameters->template estimateParameter<estimators::AbstractEstimator<double>>(estimator);
    common_intrinsics_parameters->template estimateParameter<estimators::AbstractEstimator<Eigen::VectorXd>>(estimator);
    scene::DynamicDivisionModelScene scene(cameras, stereo_pairs);
    scene.estimateStereoPairs(estimator);
    
    scene.saveScene(archivers);


    return 0;
}