#ifndef _SZ_REGRESSION_PREDICTOR_HPP
#define _SZ_REGRESSION_PREDICTOR_HPP

#include "def.hpp"
#include "utils/Iterator.hpp"
#include "predictor/Predictor.hpp"
#include "quantizer/IntegerQuantizer.hpp"
#include "encoder/HuffmanEncoder.hpp"
#include <cstring>
#include <iostream>
#include <fstream>

namespace SZ {

// N-d regression predictor
    template<class T>
    class RegressionPredictor : public concepts::PredictorInterface<T> {
    public:
        static const uint8_t predictor_id = 0b00000010;
        size_t N;
        RegressionPredictor(size_t N) : quantizer_independent(0), quantizer_liner(0), prev_coeffs(std::vector<T>(N+1)), current_coeffs(std::vector<T>(N+1)), N(N) {}

        RegressionPredictor(uint block_size, T eb, size_t N) : quantizer_independent(eb / (N + 1)),
                                                     quantizer_liner(eb / (N + 1) / block_size),
                                                     prev_coeffs(std::vector<T>(N+1)), current_coeffs(std::vector<T>(N+1)), N(N) {
        }

        RegressionPredictor(uint block_size, T eb1, T eb2, size_t N) : quantizer_independent(eb1),
                                                             quantizer_liner(eb2),
                                                             prev_coeffs(std::vector<T>(N+1)), current_coeffs(std::vector<T>(N+1)),N(N) {
        }

        using Range = multi_dimensional_range<T>;
        using iterator = typename multi_dimensional_range<T>::iterator;

        void precompress_data(const iterator &) const noexcept {}

        void postcompress_data(const iterator &) const noexcept {}

        void predecompress_data(const iterator &) const noexcept {}

        void postdecompress_data(const iterator &) const noexcept {}

        inline T estimate_error(const iterator &iter) const noexcept {
            return fabs(*iter - predict(iter));
        }

        bool precompress_block(const std::shared_ptr<Range> &range) noexcept {
            // std::cout << "precompress_block" << std::endl;
            auto dims = range->get_dimensions();
            size_t num_elements = 1;
            for (const auto &dim : dims) {
                num_elements *= dim;
                if (dim <= 1) {
                    return false;
                }
            }

            T num_elements_recip = 1.0 / num_elements;
            std::vector<double> sum(N+1);

            {
                auto range_begin = range->begin();
                auto range_end = range->end();
                for (auto iter = range_begin; iter != range_end; ++iter) {
                    double sum_cumulative = 0;
                    for (int t = 0; t < dims[N - 1]; t++) {
                        T data = *iter;
                        sum_cumulative += data;
                        sum[N - 1] += iter.get_local_index(N - 1) * data;
                        iter.move();
                    }
                    for (int i = 0; i < N - 1; i++) {
                        sum[i] += sum_cumulative * iter.get_local_index(i);
                    }
                    sum[N] += sum_cumulative;
                }
            }

            std::fill(current_coeffs.begin(), current_coeffs.end(), 0);
            current_coeffs[N] = sum[N] * num_elements_recip;
            for (int i = 0; i < N; i++) {
                current_coeffs[i] = (2 * sum[i] / (dims[i] - 1) - sum[N]) * 6 * num_elements_recip / (dims[i] + 1);
                current_coeffs[N] -= (dims[i] - 1) * current_coeffs[i] / 2;
            }
            return true;
        }

        void precompress_block_commit() noexcept {
            pred_and_quantize_coefficients();
            std::copy(current_coeffs.begin(), current_coeffs.end(), prev_coeffs.begin());
        }

        inline T predict(const iterator &iter) const noexcept {
            T pred = 0;
            for (int i = 0; i < N; i++) {
                pred += iter.get_local_index(i) * current_coeffs[i];
            }
            pred += current_coeffs[N];
            return pred;
        }

        void save(uchar *&c) const {
            std::cout << "save regression predictor" << std::endl;
            c[0] = 0b00000010;
            c += sizeof(uint8_t);
            *reinterpret_cast<size_t *>(c) = regression_coeff_quant_inds.size();
            c += sizeof(size_t);
            if (!regression_coeff_quant_inds.empty()) {
                quantizer_independent.save(c);
                quantizer_liner.save(c);
                HuffmanEncoder<int> encoder = HuffmanEncoder<int>();
                encoder.preprocess_encode(regression_coeff_quant_inds,
                                          4 * std::max(quantizer_liner.get_radius(), quantizer_independent.get_radius()));
                encoder.save(c);
                encoder.encode(regression_coeff_quant_inds, c);
                encoder.postprocess_encode();
            }
        }

        bool predecompress_block(const std::shared_ptr<Range> &range) noexcept {
            for (const auto &dim :  range->get_dimensions()) {
                if (dim <= 1) {
                    return false;
                }
            }
            pred_and_recover_coefficients();
            return true;
        }

        void load(const uchar *&c, size_t &remaining_length) {
            //TODO: adjust remaining_length
            std::cout << "load regression predictor" << std::endl;
            c += sizeof(uint8_t);
            remaining_length -= sizeof(uint8_t);

            size_t coeff_size = *reinterpret_cast<const size_t *>(c);
            c += sizeof(size_t);
            remaining_length -= sizeof(size_t);
            if (coeff_size != 0) {

                quantizer_independent.load(c, remaining_length);
                quantizer_liner.load(c, remaining_length);
                HuffmanEncoder<int> encoder = HuffmanEncoder<int>();
                encoder.load(c, remaining_length);
                regression_coeff_quant_inds = encoder.decode(c, coeff_size);
                encoder.postprocess_decode();
                remaining_length -= coeff_size * sizeof(int);
                std::fill(current_coeffs.begin(), current_coeffs.end(), 0);
                regression_coeff_index = 0;
            }
        }

        void print() const {
            std::cout << "Regression predictor, indendent term eb = " << quantizer_independent.get_eb() << "\n";
            std::cout << "Regression predictor, linear term eb = " << quantizer_liner.get_eb() << "\n";
            int count = 0;
            int ind = regression_coeff_index ? regression_coeff_index : regression_coeff_quant_inds.size();
            std::cout << "Prev coeffs: ";
            for (const auto &c:prev_coeffs) {
                std::cout << c << " ";
            }
            std::cout << "\nCurrent coeffs: ";
            for (const auto &c:current_coeffs) {
                std::cout << c << " ";
            }
            std::cout << std::endl;
        }

        void clear() {
            quantizer_liner.clear();
            quantizer_independent.clear();
            regression_coeff_quant_inds.clear();
            regression_coeff_index = 0;
            current_coeffs = {0};
            prev_coeffs = {0};
        }

    private:
        LinearQuantizer<T> quantizer_liner, quantizer_independent;
        std::vector<int> regression_coeff_quant_inds;
        size_t regression_coeff_index = 0;
        std::vector<T> current_coeffs;
        std::vector<T> prev_coeffs;

        void pred_and_quantize_coefficients() {
            for (int i = 0; i < N; i++) {
                regression_coeff_quant_inds.push_back(quantizer_liner.quantize_and_overwrite(current_coeffs[i], prev_coeffs[i]));
            }
            regression_coeff_quant_inds.push_back(
                    quantizer_independent.quantize_and_overwrite(current_coeffs[N], prev_coeffs[N]));
        }

        void pred_and_recover_coefficients() {
            for (int i = 0; i < N; i++) {
                current_coeffs[i] = quantizer_liner.recover(current_coeffs[i],
                                                            regression_coeff_quant_inds[regression_coeff_index++]);
            }
            current_coeffs[N] = quantizer_independent.recover(current_coeffs[N],
                                                              regression_coeff_quant_inds[regression_coeff_index++]);
        }
    };

}

#endif
