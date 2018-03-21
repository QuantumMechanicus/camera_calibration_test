//
// Created by danielbord on 2/14/18.
//
#include <boost/program_options.hpp>
#include <iostream>
#include <Core.h>

int main(int argc, char *argv[]) {

    std::vector<std::string> f_infosa, f_infosb, f_infosc, f_inlal, f_inlar, f_tr;
    int number_of_distortion_coefficients;
    int non_linear_iter;
    double percent_of_inliers;

    namespace po = boost::program_options;
    try {

        po::options_description desc("Automatic solver options");
        desc.add_options()
                ("help", "Print help message")
                ("a", po::value<std::vector<std::string>>(&f_infosa)->multitoken()->required(),
                 "File with data describes path to keypoints, cameras info etc (see examples)")
                ("l", po::value<std::vector<std::string>>(&f_inlal)->multitoken()->required(),
                 "File with data describes path to keypoints, cameras info etc (see examples)")
                ("r", po::value<std::vector<std::string>>(&f_inlar)->multitoken()->required(),
                 "File with data describes path to keypoints, cameras info etc (see examples)")
                ("tr", po::value<std::vector<std::string>>(&f_tr)->multitoken()->required(),
                 "File with data describes path to keypoints, cameras info etc (see examples)")

            /*("b", po::value<std::vector<std::string>>(&f_infosb)->multitoken()->required(),
             "File with data describes path to keypoints, cameras info etc (see examples)")
            ("c", po::value<std::vector<std::string>>(&f_infosc)->multitoken()->required(),
             "File with data describes path to keypoints, cameras info etc (see examples)")*/;

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
    auto number_of_pairs = f_infosa.size();
    scene::StdVector<Sophus::SE3d> a, b, c;
    scene::StdVector<scene::WorldPoints> tr;
    scene::StdVector<scene::ImagePoints> inlal, inlar;
    for (size_t k = 0; k < number_of_pairs; ++k) {
        std::cout << f_infosa[k] << std::endl;
        std::cout << f_inlal[k] << std::endl;
        Eigen::Matrix4d r;
        scene::WorldPoints ww;
        utils::loadMatrix(f_tr[k], ww);
        tr.push_back(ww);
        scene::ImagePoints inll, inlr;
        utils::loadMatrix(f_inlal[k], inll);
        std::cout << f_inlar[k] << std::endl;
        utils::loadMatrix(f_inlar[k], inlr);
        std::cout << f_infosa[k] << std::endl;
        utils::loadMatrix(f_infosa[k], r);
        Sophus::SE3d tmp;
        tmp = Sophus::SE3d::fitToSE3(r);
        inlal.push_back(inll);
        inlar.push_back(inlr);

        a.push_back(tmp);

        /*utils::loadMatrix(f_infosb[k], r);
        tmp = Sophus::SE3d::fitToSE3(r);

        b.push_back(tmp);

        utils::loadMatrix(f_infosc[k], r);
        tmp = Sophus::SE3d::fitToSE3(r);

        c.push_back(tmp);*/

    }


    for (size_t k = 0; k < number_of_pairs; ++k)
    {
        std::cout << k + 1<< "-th pair rotation angles: \n";
        std::cout << (a[k].so3().log()).norm()*180/M_PI << std::endl;
        //std::cout << (b[k].so3().log()).norm()*180/M_PI << std::endl;
        //std::cout << (c[k].so3().log()).norm()*180/M_PI << std::endl;

    }
    for (size_t k = 0; k < number_of_pairs; ++k)
    {
        std::cout << k + 1<< "-th pair translations: \n";
        std::cout << a[k].translation().normalized().transpose() << std::endl;
        //std::cout << b[k].translation().normalized().transpose() << std::endl;
        //std::cout << c[k].translation().normalized().transpose() << std::endl;
    }
    std::vector<double> meds;
    for (size_t k = 1; k < number_of_pairs; ++k) {
        std::vector<size_t> left_ind, right_ind;
        scene::StdVector<scene::ImagePoint> both;
        scene::ImagePoints left_corr = inlar[k - 1];
        scene::ImagePoints right_corr = inlal[k];

        for (size_t i = 0; i < left_corr.cols(); ++i) {
            for (size_t j = 0; j < right_corr.cols(); ++j) {
                if ((left_corr.col(i) - right_corr.col(j)).norm() < 1e-8) {
                    left_ind.push_back(i);
                    right_ind.push_back(j);
                    both.push_back(left_corr.col(i));
                }
            }
        }
        std::vector<double> ratios;

        Eigen::Matrix3d calibration_matrix;
        calibration_matrix.setIdentity();
        calibration_matrix(0, 0) = calibration_matrix(1, 1) = 0.774162;
        Eigen::VectorXd dist(2);
        dist << -0.64759,
                -0.251316;
        for (size_t i = 0; i < both.size(); ++i) {
            scene::ImagePoint leftu = inlal[k - 1].col(left_ind[i]);
            scene::ImagePoint middleu = both[i];
            scene::ImagePoint rightu = inlar[k].col(right_ind[i]);

            Sophus::SE3d leftToMid = a[k - 1];
            Sophus::SE3d midToRight = a[k];
            leftToMid.translation() = leftToMid.translation().normalized();
            midToRight.translation() = midToRight.translation().normalized();

            Eigen::Vector3d backpojected_in_left_coors = tr[k - 1].col(left_ind[i]);
            Eigen::Vector3d backpojected_in_mid_coors = leftToMid*backpojected_in_left_coors;
            Eigen::Vector3d backpojected_in_mid_coors2 = tr[k].col(right_ind[i]);

            std::cout << "Img vs world\n";
            std::cout << utils::distortion_problem::undistortion(both[i], dist).transpose() << std::endl;
            std::cout << (calibration_matrix*backpojected_in_mid_coors).hnormalized().eval().transpose() << std::endl;
            std::cout << (calibration_matrix*backpojected_in_mid_coors2).hnormalized().eval().transpose() << std::endl;

            double dd1 =  backpojected_in_mid_coors(2),
            dd2 = backpojected_in_mid_coors2(2);
            ratios.push_back(dd1/dd2);

            std::cout <<"___: "<< backpojected_in_mid_coors.transpose() << std::endl;
            std::cout << backpojected_in_mid_coors2.transpose() << "___\n"<< std::endl;

            /*Eigen::Vector3d backpojected_in_middle_coors = utils::triangulate(leftToMid,
                                                                              calibration_matrix.inverse() * leftu.homogeneous(),
                                                                              calibration_matrix.inverse() * middleu.homogeneous());

            Eigen::Vector3d backpojected_in_left_coors = leftToMid.inverse() * backpojected_in_middle_coors;

            double dd1, dd2;

            dd1 = backpojected_in_middle_coors(2);

            Eigen::Vector3d backpojected_in_right_coors = utils::triangulate(midToRight,
                                                                             calibration_matrix.inverse() * middleu.homogeneous(),
                                                                             calibration_matrix.inverse() * rightu.homogeneous());


            Eigen::Vector3d backpojected_in_middle_coors2 = midToRight.inverse() * backpojected_in_right_coors;

            dd2 = backpojected_in_middle_coors2(2);

            std::cout << backpojected_in_middle_coors.transpose() << std::endl;
            std::cout << backpojected_in_middle_coors2.transpose() << std::endl;
            ratios.push_back(dd1 / dd2);*/
        }

        std::cout << "Ratios for " << k << "th pair:\n";
        for (size_t i = 0; i < both.size(); ++i) {
            std::cout << ratios[i] << " ";
        }
        std::cout << std::endl;
        int med = static_cast<int>(both.size() / 2);
        std::nth_element(ratios.begin(), ratios.begin() + med, ratios.end());
        std::cout << "kth median: " << ratios[med] << std::endl;
        meds.push_back(ratios[med]);
    }
    a[0].translation() = a[0].translation().normalized();

    /*for (size_t k = 0; k < number_of_pairs; ++k)
    {
        for (size_t j = 0; j < tr[k].cols(); ++j)
        {
            tr[k].col(k) = tr[k].col(k)/tr[k].col(k)(2);
        }
    }*/
    double mul = meds[0];
    utils::saveMatrix(f_infosa[0]+"new",a[0].matrix(), true);
    utils::saveMatrix(f_tr[0]+"new", tr[0], true);
    Sophus::SE3d motion;
    for (size_t k = 1; k < number_of_pairs; ++k)
    {
        tr[k] *= mul;
        a[k].translation() = a[k].translation().normalized()*mul;
        motion = motion*a[k].inverse();
        std::cout << tr[k].size() << std::endl;
        for (size_t j = 0; j < tr[k].cols(); ++j)
         tr[k].col(j) = motion*tr[k].col(j);

        std::cout << tr[k].size() << std::endl;
        utils::saveMatrix(f_infosa[k]+"new",a[k].matrix(), true);
        utils::saveMatrix(f_tr[k]+"new", tr[k], true);
        if (k < number_of_pairs - 1)
            mul *= meds[k];

    }


    /*for (size_t k = 0; k < number_of_pairs; ++k)
    {
        std::cout << k + 1<< "-th pair rotation angles: \n";
        std::cout << (a[k].so3().log()).norm()*180/M_PI << std::endl;
        std::cout << (b[k].so3().log()).norm()*180/M_PI << std::endl;
        std::cout << (c[k].so3().log()).norm()*180/M_PI << std::endl;

    }
    for (size_t k = 0; k < number_of_pairs; ++k)
    {
        std::cout << k + 1<< "-th pair translations: \n";
        std::cout << a[k].translation().normalized().transpose() << std::endl;
        std::cout << b[k].translation().normalized().transpose() << std::endl;
        std::cout << c[k].translation().normalized().transpose() << std::endl;
    }*/
    return 0;
}
