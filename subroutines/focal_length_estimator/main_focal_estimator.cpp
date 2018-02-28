//
// Created by danielbord on 2/14/18.
//
#include <boost/program_options.hpp>
#include "Focal_Estimator.h"

int main(int argc, char *argv[]) {
    std::vector<std::string> f_infos;


    namespace po = boost::program_options;
    try {

        po::options_description desc("Automatic solver options");
        desc.add_options()
                ("help", "Print help message")
                ("d", po::value<std::vector<std::string>>(&f_infos)->multitoken()->required(),
                 "File with data describes path to keypoints, cameras info etc (see examples)");

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
    double w = 0;
    double h = 0;
    auto number_of_pairs = f_infos.size();

    std::shared_ptr<scene::DynamicDivisionModelStereoPair::VertexMap_t> cameras
            = std::make_shared<scene::DynamicDivisionModelStereoPair::VertexMap_t>(
                    scene::DynamicDivisionModelStereoPair::VertexMap_t());

    std::vector<scene::DynamicDivisionModelStereoPair> stereo_pairs(number_of_pairs);

    std::vector<scene_serialization::DynamicSceneArchiver> archivers(
            number_of_pairs);
    scene::FundamentalMatrices fundamental_matrices(number_of_pairs);
    for (size_t k = 0; k < number_of_pairs; ++k) {
        
        archivers[k].parse(f_infos[k]);
        stereo_pairs[k].loadScene(archivers[k], cameras);
        fundamental_matrices[k] = stereo_pairs[k].getFundamentalMatrix();
        w += (stereo_pairs[k].getLeftIntrinsicsPointer()->getWidth() +
              stereo_pairs[k].getRightIntrinsicsPointer()->getWidth()) / 2.0;

        h += (stereo_pairs[k].getLeftIntrinsicsPointer()->getHeight() +
              stereo_pairs[k].getRightIntrinsicsPointer()->getHeight()) / 2.0;
    }
    w /= number_of_pairs;
    h /= number_of_pairs;
    std::cout << (*cameras).size() << std::endl;
    estimators::SimpleFocalEstimator estimator(fundamental_matrices, static_cast<unsigned int>(w),
                                               static_cast<unsigned int>(h));
    for (size_t k = 0; k < number_of_pairs; ++k) {
        stereo_pairs[k].estimateLeftCamera(estimator);
        if (k == number_of_pairs - 1)
        {
            stereo_pairs[number_of_pairs-1].estimateRightCamera(estimator);
        }
        stereo_pairs[k].saveScene(archivers[k]);
    }


    return 0;
}