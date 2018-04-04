//
// Created by danielbord on 3/1/18.
//

#include "Global_Non_Linear_Estimator.h"

namespace non_linear_optimization {

    GlobalNonLinearEstimatorOptions::GlobalNonLinearEstimatorOptions(double quantile_to_minimize,
                                                                     double max_interval,
                                                                     double image_radius, double w, double h) :
            max_interval_(max_interval),
            quantile_to_minimize_(quantile_to_minimize),
            image_radius_(image_radius), w_(w), h_(h) {}

    void GlobalNonLinearEstimator::estimateImpl() {
        is_estimated_ = true;
        int residuals = 0;
        long number_of_distortion_coefficients = lambdas_.rows();

        ceres::Problem problem;
        double *lambda_ptr = lambdas_.data();
        double *focal_ptr = &focal_length_;
        double pp_ptr[2];
        pp_ptr[0] = ppx_;
        pp_ptr[1] = ppy_;


        problem.AddParameterBlock(lambda_ptr, static_cast<int>(number_of_distortion_coefficients));
        //problem.SetParameterBlockConstant(lambda_ptr);
        problem.AddParameterBlock(focal_ptr, 1);
        problem.AddParameterBlock(pp_ptr, 2);
        //problem.SetParameterBlockConstant(focal_ptr);
        //problem.SetParameterBlockConstant(pp_ptr);
        Eigen::Matrix3d calibration_matrix;
        calibration_matrix.setIdentity();
        calibration_matrix(0, 0) = calibration_matrix(1, 1) = *focal_ptr;
        calibration_matrix(0, 2) = pp_ptr[0];
        calibration_matrix(1, 2) = pp_ptr[1];
        std::cout << *focal_ptr << " " << focal_length_ << " Check" << std::endl;
        scene::StdVector<scene::WorldPoints> wp(number_of_pairs_);
        scene::StdVector<scene::ImagePoints> i1d(number_of_pairs_), i2d(number_of_pairs_);
        for (size_t kth_pair = 0; kth_pair < number_of_pairs_; ++kth_pair) {
            auto &kth_translation = translations_[kth_pair];
            auto &kth_rotation = rotations_[kth_pair];
            auto &kth_fm = fundamental_matrices_[kth_pair];
            //TODO to function essintial to fund
            kth_fm = calibration_matrix.transpose().inverse() * utils::screw_hat(kth_translation) *
                     kth_rotation.matrix() * calibration_matrix.inverse();

            std::cout << kth_fm << std::endl;

            double *tr_ptr = kth_translation.data();
            double *rot_ptr = kth_rotation.data();
            problem.AddParameterBlock(tr_ptr, 3, new local_parametrization::LocalParameterizationSphere(1));
            problem.AddParameterBlock(rot_ptr, 4, new local_parametrization::LocalParameterizationSO3());


            std::vector<size_t> inliers_ind;
            Sophus::SE3d leftToRight(kth_rotation, kth_translation);
            double interval = utils::distortion_problem::findInliers(left_pictures_keypoints_[kth_pair],
                                                                     right_pictures_keypoints_[kth_pair],
                                                                     lambdas_,
                                                                     leftToRight, calibration_matrix,
                                                                     options_.quantile_to_minimize_,
                                                                     inliers_ind, options_.image_radius_);
            std::cout << "Non line Interval:: " << interval << std::endl;


            i1d[kth_pair].resize(Eigen::NoChange, inliers_ind.size());
            i2d[kth_pair].resize(Eigen::NoChange, inliers_ind.size());
            for (size_t kth_inlier = 0; kth_inlier < inliers_ind.size(); ++kth_inlier) {
                i1d[kth_pair].col(kth_inlier) = left_pictures_keypoints_[kth_pair].col(inliers_ind[kth_inlier]);
                i2d[kth_pair].col(kth_inlier) = right_pictures_keypoints_[kth_pair].col(inliers_ind[kth_inlier]);
            }
            wp[kth_pair].resize(Eigen::NoChange, i1d[kth_pair].cols());

            LOG(INFO) << kth_pair << " : " << i1d[kth_pair].cols() << " inliers [out of "
                      << left_pictures_keypoints_[kth_pair].cols() << "]";
            for (size_t k = 0; k < i1d[kth_pair].cols(); ++k) {
                auto kth_translation_matrix = utils::screw_hat(kth_translation);
                scene::WorldPoint right_bprj;
                scene::WorldPoint left_bprj;
                std::cout << "Test:: " << calibration_matrix.transpose() * kth_fm * calibration_matrix -
                                          kth_translation_matrix * kth_rotation.matrix() << std::endl << std::endl;
                scene::ImagePoint left_c, right_c;

                kth_fm = calibration_matrix.inverse().transpose() * kth_translation_matrix * kth_rotation.matrix() *
                         calibration_matrix.inverse();

                std::cout << "Calib:: " << calibration_matrix << std::endl;
                std::cout << "Lamb: " << lambdas_.transpose() << std::endl;
                std::cout << "Dos " << utils::distortion_problem::undistortion<double>(i1d[kth_pair].col(k),
                                                                                       lambdas_) << std::endl;
                utils::triangulate(kth_fm, kth_rotation.matrix(), kth_translation_matrix,
                                   utils::distortion_problem::undistortion<double>(i1d[kth_pair].col(k),
                                                                                   lambdas_),
                                   utils::distortion_problem::undistortion<double>(i2d[kth_pair].col(k), lambdas_),
                                   left_bprj, right_bprj, lambdas_.transpose(), calibration_matrix);
                wp[kth_pair].col(k) = left_bprj;


                std::cout << (utils::distortion_problem::distortion<double>(
                        (calibration_matrix * right_bprj).hnormalized().eval(), lambdas_)
                              - i2d[kth_pair].col(k)).transpose() * options_.image_radius_ << std::endl;

                std::cout << (utils::distortion_problem::distortion<double>(
                        (calibration_matrix * left_bprj).hnormalized().eval(), lambdas_)
                              - i1d[kth_pair].col(k)).transpose() * options_.image_radius_ << std::endl;

                double *wp_ptr = &wp[kth_pair](0, k);

                problem.AddParameterBlock(wp_ptr, 3);


                auto fun = new ceres::DynamicAutoDiffCostFunction<GlobalOptimizerFunctor<>>(
                        new GlobalOptimizerFunctor<>(i1d[kth_pair].col(k), i2d[kth_pair].col(k),
                                                     static_cast<int>(number_of_distortion_coefficients),
                                                     options_.image_radius_));
                fun->AddParameterBlock(static_cast<int>(number_of_distortion_coefficients));
                fun->AddParameterBlock(1);
                fun->AddParameterBlock(2);
                fun->AddParameterBlock(3);
                fun->AddParameterBlock(4);
                fun->AddParameterBlock(3);
                fun->SetNumResiduals(4);
                problem.AddResidualBlock(fun, new ceres::HuberLoss(5.0), lambda_ptr, focal_ptr, pp_ptr,
                                         tr_ptr, rot_ptr, wp_ptr);

                ++residuals;

                double *inputs[] = {lambda_ptr, focal_ptr, pp_ptr, tr_ptr, rot_ptr, wp_ptr};
                Eigen::Vector4d error;
                fun->Evaluate(inputs, error.data(), nullptr);
                std::cout << *wp_ptr << std::endl;
                LOG(INFO) << "Functor returns error: " << error.transpose();
            }

        }
        LOG(INFO) << residuals << " points";
        LOG(INFO) << lambdas_ << " --- coefficients before estimation";
        ceres::Solver::Options options;
        options.max_num_iterations = 500;
        options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
        options.num_threads = 8;
        options.num_linear_solver_threads = 8;
        options.function_tolerance = 1e-16;
        options.parameter_tolerance = 1e-16;
        options.minimizer_progress_to_stdout = true;


        // Solve
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);
        LOG(INFO) << summary.BriefReport();
        LOG(INFO) << lambdas_ << " --- estimated coefficients";
        LOG(INFO) << focal_length_ << " --- focal length";
        LOG(INFO) << ppx_ << " " << ppy_;
        LOG(INFO) << 2 * atan(1 / focal_length_) * 180 / M_PI;
        LOG(INFO) << std::sqrt(summary.final_cost * 2 / residuals) << " final rmse";

        const int FOVS = 6;
        double w = options_.w_;
        double h = options_.h_;
        double r = options_.image_radius_;
        scene::ImagePoint center(w / 2, h / 2);
        double fov[FOVS];
        scene::ImagePoint points[FOVS][2] = {
                {Eigen::Vector2d(0, 0),       Eigen::Vector2d(w, h)},
                {Eigen::Vector2d(0, h / 2),   Eigen::Vector2d(w, h / 2)},
                {Eigen::Vector2d(w / 2, 0),   Eigen::Vector2d(w / 2, h)},
                {Eigen::Vector2d(444, h / 2), Eigen::Vector2d(6486, h / 2)},
                {Eigen::Vector2d(35, h / 2),  Eigen::Vector2d(6672, h / 2)},
                {Eigen::Vector2d(176, h / 2), Eigen::Vector2d(6813, h / 2)}
        };

        std::string tag[FOVS] = {
                "DFOV",
                "HFOV",
                "VFOV",
                "MAGICA",
                "MAGICB",
                "MAGICC"
        };

        std::cout << "Cal mat " << calibration_matrix << std::endl;
        calibration_matrix(0, 0) = calibration_matrix(1, 1) = *focal_ptr;

        for (int i = 0; i < FOVS; ++i) {
            auto undistorted1 = utils::distortion_problem::undistortion<double>((points[i][0] - center) / r, lambdas_);
            auto undistorted2 = utils::distortion_problem::undistortion<double>((points[i][1] - center) / r, lambdas_);

            Eigen::Vector3d ray1 = (calibration_matrix.inverse() * undistorted1.homogeneous()).normalized(),
                    ray2 = (calibration_matrix.inverse() * undistorted2.homogeneous()).normalized();
            LOG(INFO) << tag[i] << ": " << (fov[i] = std::acos(ray1.dot(ray2)) * 180 / M_PI);


        }
        double total = 0.0;
        for (int i = FOVS - 3; i < FOVS; ++i)
            total += fov[i];
        LOG(INFO) << "Total horizontal FOV guesstimate: " << total << " degrees";
        for (size_t kth_pair = 0; kth_pair < number_of_pairs_; ++kth_pair) {
            calibration_matrix.setIdentity();
            calibration_matrix(0, 0) = calibration_matrix(1, 1) = *focal_ptr;
            calibration_matrix(0, 2) = pp_ptr[0];
            calibration_matrix(1, 2) = pp_ptr[1];
            fundamental_matrices_[kth_pair] =
                    calibration_matrix.inverse().transpose() * utils::screw_hat(translations_[kth_pair]) *
                    rotations_[kth_pair].matrix() * calibration_matrix.inverse();
            auto &kth_fundamental_matrix = fundamental_matrices_[kth_pair];
            Eigen::JacobiSVD<Eigen::Matrix3d> fmatrix_svd(kth_fundamental_matrix,
                                                          Eigen::ComputeFullU | Eigen::ComputeFullV);
            Eigen::Vector3d singular_values = fmatrix_svd.singularValues();
            singular_values[2] = 0.0;
            kth_fundamental_matrix = fmatrix_svd.matrixU() * singular_values.asDiagonal() *
                                     fmatrix_svd.matrixV().transpose();
            kth_fundamental_matrix /= kth_fundamental_matrix(2, 2);
            
            std::cout << "Save: " << kth_pair << std::endl;
            utils::saveMatrix("./left_inl_pair_" + std::to_string(kth_pair), i1d[kth_pair], true);
            utils::saveMatrix("./right_inl_pair_" + std::to_string(kth_pair), i2d[kth_pair], true);
            scene::ImagePoints i1un, i2un;
            i1un = i1d[kth_pair];
            i2un = i2d[kth_pair];
            std::cout << "Lambdas: " << lambdas_.transpose() << std::endl;
            double d = std::max(w, h) / 2.0;
            double dr = d / r;
            double alpha = utils::distortion_problem::undistortionDenominator<double>(dr, lambdas_.cast<double>());

            for (size_t ll = 0; ll < i1d[kth_pair].cols(); ++ll) {
                i1un.col(ll) = alpha*r*utils::distortion_problem::undistortion<double>(i1d[kth_pair].col(ll), lambdas_) + center;
                i2un.col(ll) = alpha*r*utils::distortion_problem::undistortion<double>(i2d[kth_pair].col(ll), lambdas_) + center;

            }
            Eigen::Matrix3d recF = fundamental_matrices_[kth_pair];
            Eigen::Matrix3d remove_center, remove_scale;
            remove_center << 1, 0, -w/2,
                             0, 1, -h/2,
                             0, 0, 1;
            remove_scale << 1.0/(alpha*r), 0, 0,
                            0, 1.0/(alpha*r), 0,
                            0, 0, 1;
            recF = remove_center.transpose()*remove_scale.transpose()*recF*remove_scale*remove_center;
            utils::saveMatrix("./left_inl_pair_und" + std::to_string(kth_pair), i1un, true);
            utils::saveMatrix("./right_inl_pair_und" + std::to_string(kth_pair), i2un, true);
            utils::saveMatrix("./triang_" + std::to_string(kth_pair), wp[kth_pair], true);
            utils::saveMatrix("./rec_F_" + std::to_string(kth_pair), recF, true);
            //LOG(INFO) << "Triangulated points optimized, " << kth_pair <<"-th\n\n" << wp[kth_pair] << std::endl << std::endl;

        }


    }


    void GlobalNonLinearEstimator::getEstimationImpl(Eigen::VectorXd &result) {
        result = lambdas_;
    }

    void GlobalNonLinearEstimator::getEstimationImpl(intrinsics::FocalLength &result) {
        result = focal_length_;
    }

    void GlobalNonLinearEstimator::getEstimationImpl(scene::FundamentalMatrices &result) {
        for (auto &kth_fundamental_matrix : fundamental_matrices_) {
            Eigen::JacobiSVD<Eigen::Matrix3d> fmatrix_svd(kth_fundamental_matrix,
                                                          Eigen::ComputeFullU | Eigen::ComputeFullV);
            Eigen::Vector3d singular_values = fmatrix_svd.singularValues();
            singular_values[2] = 0.0;
            kth_fundamental_matrix = fmatrix_svd.matrixU() * singular_values.asDiagonal() *
                                     fmatrix_svd.matrixV().transpose();
            kth_fundamental_matrix /= kth_fundamental_matrix(2, 2);

        }
        result = fundamental_matrices_;
    }


    bool GlobalNonLinearEstimator::isEstimated() const {
        return is_estimated_;
    }

    GlobalNonLinearEstimator::GlobalNonLinearEstimator(
            const scene::StdVector<scene::ImagePoints> &left_pictures_keypoints,
            const scene::StdVector<scene::ImagePoints> &right_pictures_keypoints,
            const Eigen::VectorXd &lambdas,
            const scene::StdVector<Sophus::SO3d> &rotations,
            const scene::StdVector<Eigen::Vector3d> &translations,
            double focal_length, double ppx, double ppy,
            GlobalNonLinearEstimatorOptions options)
            : rotations_((rotations)), translations_((translations)),
              focal_length_(focal_length), ppx_(ppx), ppy_(ppy),
              lambdas_((lambdas)), left_pictures_keypoints_((left_pictures_keypoints)),
              right_pictures_keypoints_((right_pictures_keypoints)), options_(options),
              is_estimated_(false) {
        number_of_pairs_ = rotations_.size();
        fundamental_matrices_.resize(number_of_pairs_);

    }


}