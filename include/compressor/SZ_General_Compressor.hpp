#ifndef _SZ_SZ_GENERAL_HPP
#define _SZ_SZ_GENERAL_HPP

#include "predictor/Predictor.hpp"
#include "predictor/LorenzoPredictor.hpp"
#include "quantizer/Quantizer.hpp"
#include "encoder/Encoder.hpp"
#include "lossless/Lossless.hpp"
#include "utils/Iterator.hpp"
#include "utils/MemoryOps.hpp"
#include "utils/Config.hpp"
#include "utils/FileUtil.h"
#include "def.hpp"
#include <cstring>
#include <encoder/HuffmanEncoder.hpp>

namespace SZ {
    template<class T, class Predictor, class Quantizer, class Encoder, class Lossless>
    class SZ_General_Compressor {
    public:
        SZ_General_Compressor(const Config<T> &conf,
                              Predictor predictor, Quantizer quantizer, Encoder encoder, Lossless lossless) :
                fallback_predictor(LorenzoPredictor<T,1>(conf.eb, conf.N)),
                predictor(predictor), quantizer(quantizer), encoder(encoder), lossless(lossless),
                block_size(conf.block_size), stride(conf.stride),
                global_dimensions(conf.dims), num_elements(conf.num),
                my_lorenzo(LorenzoPredictor<T,1>(conf.eb, conf.N)),N(conf.N){
            prediction_debug = std::vector<T>(num_elements);
            element_debug = std::vector<T>(num_elements);
            static_assert(std::is_base_of_v<concepts::PredictorInterface<T>, Predictor>,
                          "must implement the predictor interface");
            static_assert(std::is_base_of_v<concepts::QuantizerInterface<T>, Quantizer>, "must implement the quatizer interface");
            static_assert(std::is_base_of_v<concepts::EncoderInterface<int>, Encoder>, "must implement the encoder interface");
            static_assert(std::is_base_of_v<concepts::LosslessInterface, Lossless>, "must implement the lossless interface");
        }

        uchar *compress_withBG(T *data, size_t &compressed_size, T bg_data, T low_range, T high_range, bool use_bitmap=false, bool preserve_sign=false, bool has_bg=false) {
            auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<T>>(data,
                                                                                         std::begin(global_dimensions),
                                                                                         std::end(global_dimensions),
                                                                                         stride, 0, N);
            auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<T>>(data,
                                                                                         std::begin(global_dimensions),
                                                                                         std::end(global_dimensions), 1,
                                                                                         0, N);

            std::vector<int> quant_inds(num_elements);
            std::vector<int> bitmap;
            std::vector<int> signs;
            if(has_bg) {
                if (preserve_sign) {
                    signs.resize(num_elements);
                }
                if (use_bitmap) {
                    bitmap.resize(num_elements);
                }

                T tmp;
                for (auto element = intra_block_range->begin(); element != intra_block_range->end(); ++element) {
                    if (*element == bg_data) {
                        if (use_bitmap) {
                            bitmap[element.get_offset()] = 1;
                        } else {
                            quant_inds[element.get_offset()] = 2 * quantizer.get_radius() + 1;
                        }
                        tmp = my_lorenzo.predict(element);
                        if (tmp >= low_range && tmp <= high_range) {
                            *element = tmp;
                        } else {
//                        tmp = (low_range + high_range)/2;
                            tmp = 0;
                            *element = tmp;
                        }
                    }
                    if (preserve_sign) {
                        if (*element > 0) {
                            signs[element.get_offset()] = 1;
                        } else if (*element == 0) {
                            signs[element.get_offset()] = 0;
                        } else if (*element < 0) {
                            signs[element.get_offset()] = 2;
                        }
                    }
                }
            }
//            if(use_bitmap) {
//                convertIntArray2ByteArray_fast_1b(bitmap, num_elements, result);
//                delete[] bitmap;
//            }

//            SZ::writefile("/Users/apple/Development/globus/cleanBG.bin", (float*)data, num_elements);
//            printf("File cleanBG saved!\n");
            std::vector<size_t> intra_block_dims(N);
            predictor.precompress_data(inter_block_range->begin());
            quantizer.precompress_data();
            size_t quant_count = 0;
            struct timespec start, end;
            clock_gettime(CLOCK_REALTIME, &start);
            {
                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();
                for (auto block = inter_begin; block != inter_end; ++block) {

                    // std::cout << *block << " " << lp.predict(block) << std::endl;
                    for (int i = 0; i < intra_block_dims.size(); i++) {
                        size_t cur_index = block.get_local_index(i);
                        size_t dims = inter_block_range->get_dimensions(i);
                        intra_block_dims[i] = (cur_index == dims - 1 &&
                                               global_dimensions[i] - cur_index * stride < block_size) ?
                                              global_dimensions[i] - cur_index * stride : block_size;
                    }

                    intra_block_range->set_dimensions(intra_block_dims.begin(), intra_block_dims.end());
                    intra_block_range->set_offsets(block.get_offset());
                    intra_block_range->set_starting_position(block.get_local_index());
                    concepts::PredictorInterface<T> *predictor_withfallback = &predictor;
                    if (!predictor.precompress_block(intra_block_range)) {
                        predictor_withfallback = &fallback_predictor;
                    }
                    predictor_withfallback->precompress_block_commit();
//                    quantizer.precompress_block();
                    auto intra_begin = intra_block_range->begin();
                    auto intra_end = intra_block_range->end();
                    for (auto element = intra_begin; element != intra_end; ++element) {
                        int offset = element.get_offset();
                        if (has_bg && !use_bitmap && quant_inds[offset] == 2 * quantizer.get_radius() + 1) {
                            *element = predictor_withfallback->predict(element);
                        } else {
                            T pred = predictor_withfallback->predict(element);
                            prediction_debug[offset] = pred;
                            quant_inds[offset] = quantizer.quantize_and_overwrite(
                                    *element, pred);
                            element_debug[offset] = *element;
                            quant_count++;
                        }
                    }
                }
            }

            clock_gettime(CLOCK_REALTIME, &end);
            std::cout << "Predition & Quantization time = "
                      << (double) (end.tv_sec - start.tv_sec) +
                         (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000
                      << "s" << std::endl;

//            int z =1;
//            for (int y=0; y< 100; y++) {
//                for(int x=0;x<500;x++){
//                    printf("%g  ",data[x + 500*(y+500*z)]);
//                }
//                printf("\n");
//            }

            predictor.postcompress_data(inter_block_range->begin());
            quantizer.postcompress_data();

            uchar *compressed_data = new uchar[2 * num_elements * sizeof(T)];
            uchar *compressed_data_pos = compressed_data;
            write(global_dimensions.data(), N, compressed_data_pos);
            write(block_size, compressed_data_pos);
            predictor.save(compressed_data_pos);
            quantizer.save(compressed_data_pos);

            if(has_bg && use_bitmap) {
                HuffmanEncoder<int> bitmap_encoder = HuffmanEncoder<int>();
                bitmap_encoder.preprocess_encode(bitmap,
                                                 2);
                bitmap_encoder.save(compressed_data_pos);
                bitmap_encoder.encode(bitmap, compressed_data_pos);
                bitmap_encoder.postprocess_encode();
            }

            if(preserve_sign){
                HuffmanEncoder<int> signs_encoder = HuffmanEncoder<int>();
                signs_encoder.preprocess_encode(signs,
                                                 3);
                signs_encoder.save(compressed_data_pos);
                signs_encoder.encode(signs, compressed_data_pos);
                signs_encoder.postprocess_encode();
            }

            encoder.preprocess_encode(quant_inds, 2 * quantizer.get_radius()+2);
            encoder.save(compressed_data_pos);
            encoder.encode(quant_inds, compressed_data_pos);
            encoder.postprocess_encode();

            uchar *lossless_data = lossless.compress(compressed_data,
                                                     compressed_data_pos - compressed_data,
                                                     compressed_size);
            lossless.postcompress_data(compressed_data);
            return lossless_data;
        }

        T *decompress_withBG(uchar const *lossless_compressed_data, const size_t length, T bg, T low_range, T high_range, bool use_bitmap=false, bool preserve_sign=false,bool has_bg=false) {
            size_t remaining_length = length;
            auto compressed_data = lossless.decompress(lossless_compressed_data, remaining_length);
            uchar const *compressed_data_pos = compressed_data;
            read(global_dimensions.data(), N, compressed_data_pos, remaining_length);
            num_elements = 1;
            for (const auto &d : global_dimensions) {
                num_elements *= d;
                std::cout << d << " ";
            }
            std::cout << std::endl;
            read(block_size, compressed_data_pos, remaining_length);
            stride = block_size;
            predictor.load(compressed_data_pos, remaining_length);
            quantizer.load(compressed_data_pos, remaining_length);
            std::vector<int> bitmap;
            std::vector<int> signs;
            if(has_bg && use_bitmap) {
                HuffmanEncoder<int> bitmap_encoder = HuffmanEncoder<int>();
                bitmap_encoder.load(compressed_data_pos, remaining_length);
                bitmap = bitmap_encoder.decode(compressed_data_pos, num_elements);
            }
            if(preserve_sign){
                HuffmanEncoder<int> signs_encoder = HuffmanEncoder<int>();
                signs_encoder.load(compressed_data_pos, remaining_length);
                signs = signs_encoder.decode(compressed_data_pos, num_elements);
            }
            encoder.load(compressed_data_pos, remaining_length);

            auto quant_inds = encoder.decode(compressed_data_pos, num_elements);
            encoder.postprocess_decode();
            lossless.postdecompress_data(compressed_data);

            int const *quant_inds_pos = (int const *) quant_inds.data();
            std::vector<size_t> intra_block_dims(N);
            auto dec_data = std::make_unique<T[]>(num_elements);
            auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<T>>(dec_data.get(),
                                                                                         std::begin(global_dimensions),
                                                                                         std::end(global_dimensions),
                                                                                         block_size,
                                                                                         0,N);

            auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<T>>(dec_data.get(),
                                                                                         std::begin(global_dimensions),
                                                                                         std::end(global_dimensions), 1,
                                                                                         0,N);

            predictor.predecompress_data(inter_block_range->begin());
            quantizer.predecompress_data();
            T tmp;
            std::cout << "start decompression" << std::endl;
            {
                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();
                for (auto block = inter_begin; block != inter_end; block++) {
                    for (int i = 0; i < intra_block_dims.size(); i++) {
                        size_t cur_index = block.get_local_index(i);
                        size_t dims = inter_block_range->get_dimensions(i);
                        intra_block_dims[i] = (cur_index == dims - 1) ? global_dimensions[i] - cur_index * block_size
                                                                      : block_size;
                    }
                    intra_block_range->set_dimensions(intra_block_dims.begin(), intra_block_dims.end());
                    intra_block_range->set_offsets(block.get_offset());
                    intra_block_range->set_starting_position(block.get_local_index());

                    concepts::PredictorInterface<T> *predictor_withfallback = &predictor;
                    if (!predictor.predecompress_block(intra_block_range)) {
                        predictor_withfallback = &fallback_predictor;
                    }
                    auto intra_begin = intra_block_range->begin();
                    auto intra_end = intra_block_range->end();
                    for (auto element = intra_begin; element != intra_end; ++element) {
                        int offset = element.get_offset();
                        if(!use_bitmap && quant_inds[offset]==2 * quantizer.get_radius() +1){
                            *element = predictor_withfallback->predict(element);
//                            *element = bg;
                        }else {
                            T pred = predictor_withfallback->predict(element);
                            if(fabs(pred - prediction_debug[offset]) > 0.02){
                                // std::cerr << "Prediction inconsistent! Offset:" << offset << "Pred="<<pred <<"debug_pred="<<prediction_debug[offset]<<std::endl;
                                printf("offset: %d, pred=%.4f, debug_pred=%.4f\n", offset, pred, prediction_debug[offset]);
                                exit(1);
                            }
                            *element = quantizer.recover(pred,
                                                         quant_inds[element.get_offset()]);
                            if( fabs(*element - element_debug[offset]) > 0.02) {
                                printf("offset: %d, *element=%.4f, debug_element=%.4f\n", offset, *element, element_debug[offset]);
                                printf("offset: %d, pred=%.4f, debug_pred=%.4f\n", offset, pred, prediction_debug[offset]);
                                exit(2);
                            }
                        }
                    }
                }
            }
            if(has_bg) {
                auto intra_range = std::make_shared<SZ::multi_dimensional_range<T>>(dec_data.get(),
                                                                                       std::begin(global_dimensions),
                                                                                       std::end(global_dimensions), 1,
                                                                                       0,N);
                for (auto element = intra_range->begin(); element != intra_range->end(); element++) {
                    int index = element.get_offset();
                    bool isBG = false;
                    if (use_bitmap) {
                        if (bitmap[index] == 1) {
                            *element = bg;
                            isBG = true;
                        }
                    } else if (quant_inds[index] == 2 * quantizer.get_radius() + 1) {
                        *element = bg;
                        isBG = true;
                    }
                    if (preserve_sign && !isBG) {
                        if (signs[index] == 0) {
                            *element = 0;
                        }
                        if ((*element > 0 && signs[index] == 2) || (*element < 0 && signs[index == 1])) {
                            *element = -*element;
                        }
                    }
                }
            }
            predictor.postdecompress_data(inter_block_range->begin());
            quantizer.postdecompress_data();
            return dec_data.release();
        }



        uchar *compress(T *data, size_t &compressed_size) {

            auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<T>>(data,
                                                                                         std::begin(global_dimensions),
                                                                                         std::end(global_dimensions), stride, 0,N);
            auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<T>>(data,
                                                                                         std::begin(global_dimensions),
                                                                                         std::end(global_dimensions), 1, 0,N);
            std::vector<size_t> intra_block_dims;
            std::vector<int> quant_inds(num_elements);
            predictor.precompress_data(inter_block_range->begin());
            quantizer.precompress_data();
            size_t quant_count = 0;
            struct timespec start, end;
            clock_gettime(CLOCK_REALTIME, &start);
            {
                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();
                for (auto block = inter_begin; block != inter_end; ++block) {

                    // std::cout << *block << " " << lp.predict(block) << std::endl;
                    for (int i = 0; i < intra_block_dims.size(); i++) {
                        size_t cur_index = block.get_local_index(i);
                        size_t dims = inter_block_range->get_dimensions(i);
                        intra_block_dims[i] = (cur_index == dims - 1 && global_dimensions[i] - cur_index * stride < block_size) ?
                                              global_dimensions[i] - cur_index * stride : block_size;
                    }

                    intra_block_range->set_dimensions(intra_block_dims.begin(), intra_block_dims.end());
                    intra_block_range->set_offsets(block.get_offset());
                    intra_block_range->set_starting_position(block.get_local_index());
                    concepts::PredictorInterface<T> *predictor_withfallback = &predictor;
                    if (!predictor.precompress_block(intra_block_range)) {
                        predictor_withfallback = &fallback_predictor;
                    }
                    predictor_withfallback->precompress_block_commit();
//                    quantizer.precompress_block();
                    auto intra_begin = intra_block_range->begin();
                    auto intra_end = intra_block_range->end();
                    for (auto element = intra_begin; element != intra_end; ++element) {
                        quant_inds[quant_count++] = quantizer.quantize_and_overwrite(
                                *element, predictor_withfallback->predict(element));
                    }
                }
            }

            clock_gettime(CLOCK_REALTIME, &end);
            std::cout << "Predition & Quantization time = "
                      << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000
                      << "s" << std::endl;

            predictor.postcompress_data(inter_block_range->begin());
            quantizer.postcompress_data();

            uchar *compressed_data = new uchar[2 * num_elements * sizeof(T)];
            uchar *compressed_data_pos = compressed_data;

            write(global_dimensions.data(), N, compressed_data_pos);
            write(block_size, compressed_data_pos);
            predictor.save(compressed_data_pos);
            quantizer.save(compressed_data_pos);

            encoder.preprocess_encode(quant_inds, 4 * quantizer.get_radius());
            encoder.save(compressed_data_pos);
            encoder.encode(quant_inds, compressed_data_pos);
            encoder.postprocess_encode();

            uchar *lossless_data = lossless.compress(compressed_data,
                                                     compressed_data_pos - compressed_data,
                                                     compressed_size);
            lossless.postcompress_data(compressed_data);
            return lossless_data;
        }


        T *decompress(uchar const *lossless_compressed_data, const size_t length) {
            auto compressed_data = lossless.decompress(lossless_compressed_data, length);
            uchar const *compressed_data_pos = compressed_data;
            size_t remaining_length = length;
            read(global_dimensions.data(), N, compressed_data_pos, remaining_length);
            num_elements = 1;
            for (const auto &d : global_dimensions) {
                num_elements *= d;
                std::cout << d << " ";
            }
            std::cout << std::endl;
            read(block_size, compressed_data_pos, remaining_length);
            stride = block_size;
            predictor.load(compressed_data_pos, remaining_length);
            quantizer.load(compressed_data_pos, remaining_length);
            encoder.load(compressed_data_pos, remaining_length);

            auto quant_inds = encoder.decode(compressed_data_pos, num_elements);
            encoder.postprocess_decode();
            lossless.postdecompress_data(compressed_data);

            int const *quant_inds_pos = (int const *) quant_inds.data();
            std::vector<size_t> intra_block_dims(N);
            auto dec_data = std::make_unique<T[]>(num_elements);
            auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<T>>(dec_data.get(),
                                                                                         std::begin(global_dimensions),
                                                                                         std::end(global_dimensions), block_size,
                                                                                         0,N);

            auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<T>>(dec_data.get(),
                                                                                         std::begin(global_dimensions),
                                                                                         std::end(global_dimensions), 1, 0,N);

            predictor.predecompress_data(inter_block_range->begin());
            quantizer.predecompress_data();

            std::cout << "start decompression" << std::endl;
            {
                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();
                for (auto block = inter_begin; block != inter_end; block++) {
                    for (int i = 0; i < intra_block_dims.size(); i++) {
                        size_t cur_index = block.get_local_index(i);
                        size_t dims = inter_block_range->get_dimensions(i);
                        intra_block_dims[i] = (cur_index == dims - 1) ? global_dimensions[i] - cur_index * block_size
                                                                      : block_size;
                    }
                    intra_block_range->set_dimensions(intra_block_dims.begin(), intra_block_dims.end());
                    intra_block_range->set_offsets(block.get_offset());
                    intra_block_range->set_starting_position(block.get_local_index());

                    concepts::PredictorInterface<T> *predictor_withfallback = &predictor;
                    if (!predictor.predecompress_block(intra_block_range)) {
                        predictor_withfallback = &fallback_predictor;
                    }
                    auto intra_begin = intra_block_range->begin();
                    auto intra_end = intra_block_range->end();
                    for (auto element = intra_begin; element != intra_end; ++element) {
                        *element = quantizer.recover(predictor_withfallback->predict(element), *(quant_inds_pos++));
                    }
                }
            }
            predictor.postdecompress_data(inter_block_range->begin());
            quantizer.postdecompress_data();
            return dec_data.release();
        }

        size_t getNumElements(){
            return num_elements;
        }


    private:
        Predictor predictor;
        LorenzoPredictor<T, 1> my_lorenzo;
        LorenzoPredictor<T, 1> fallback_predictor;
        Quantizer quantizer;
        Encoder encoder;
        Lossless lossless;
        uint block_size;
        uint stride;
        size_t num_elements;
        std::vector<size_t> global_dimensions;
        std::vector<T> prediction_debug, element_debug;
        size_t N;
    };


    template<class T, class Predictor, class Quantizer, class Encoder, class Lossless>
    SZ_General_Compressor<T, Predictor, Quantizer, Encoder, Lossless>
    make_sz_general_compressor(const Config<T> &conf, Predictor predictor, Quantizer quantizer, Encoder encoder,
                               Lossless lossless, size_t N) {
        return SZ_General_Compressor<T, Predictor, Quantizer, Encoder, Lossless>(conf, predictor, quantizer, encoder,
                                                                                    lossless, N);
    }



}
#endif

