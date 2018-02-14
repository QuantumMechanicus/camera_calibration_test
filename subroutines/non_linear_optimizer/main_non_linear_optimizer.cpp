//
// Created by danielbord on 2/14/18.
//
#include <boost/program_options.hpp>
#include "Core.h"

int main(int argc, char *argv[]) {

    unsigned int w, h;
    std::vector<std::string> input1, input2, f_estimated_fundamental_matrices;
    int non_linear_iter, number_of_distortion_coefficients;
    double percent_of_inliers;
    std::string f_optimizer_results;
    namespace po = boost::program_options;
    try {

        po::options_description desc("Optimization input options");
        desc.add_options()
                ("help", "Print help message")
                ("lp", po::value<std::vector<std::string> >(&input1)->multitoken()->required(),
                 "Path to files with first (left) camera keypoints")
                ("rp", po::value<std::vector<std::string> >(&input2)->multitoken()->required(),
                 "Path to files with second (right) camera keypoints")
                ("ff",
                 po::value<std::vector<std::string> >(&f_estimated_fundamental_matrices)->multitoken()->required(),
                 "File with %n_pic estimated fundamental matrices")
                ("rf", po::value<std::string>(&f_optimizer_results)->default_value("./OptimizerResults/"),
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

    assert(input1.size() == input2.size() && "number of left and right cameras should be equal");
    assert(f_estimated_fundamental_matrices.size() == input1.size() &&
           "number of cameras and fundamental matrices should be equal");



    return 0;
}