//
// Created by Kai Zhao on 4/28/20.
//

#ifndef SZ_CONFIG_HPP
#define SZ_CONFIG_HPP

namespace SZ {
    template<class T>
    class Config {
    public:
        Config(T _eb, std::vector<size_t> _dims) : eb(_eb), dims(_dims) {
            switch (_dims.size()) {
                case 1:
                    block_size = 128;
                    break;
                case 2:
                    block_size = 16;
                    break;
                default:
                    // >= 3D
                    block_size = 6;
                    break;
            }
            stride = block_size;
            num = 1;
            N=_dims.size();
            for (int i=0;i<N;i++) {
                num *= _dims[i];
            }
        }

        std::vector<size_t> dims;
        size_t num;
        size_t N;
        bool enable_lorenzo = true;
        bool enable_2ndlorenzo = false;
        bool enable_regression = true;
        bool enable_2ndregression = false;
        bool enable_lossless = true;
        size_t quant_bin = 32768;
        uint block_size, stride, pred_dim = 0;
        T eb;
    };
}

#endif //SZ_CONFIG_HPP
