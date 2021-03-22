#ifndef _SZ_INTEGER_QUANTIZER_HPP
#define _SZ_INTEGER_QUANTIZER_HPP

#include <cstring>
#include <cassert>
#include <iostream>
#include <vector>
#include "def.hpp"
#include "quantizer/Quantizer.hpp"

namespace SZ {

// data with T type
// return int
    template<class T>
    class PredictionBasedQuantizer : public concepts::QuantizerInterface<T> {
    protected:
        T error_bound;
        T error_bound_reciprocal;
        int radius; // quantization interval radius
    public:
        ~PredictionBasedQuantizer() = default;

        PredictionBasedQuantizer() = default;

        PredictionBasedQuantizer(PredictionBasedQuantizer const &) = default;

        PredictionBasedQuantizer(PredictionBasedQuantizer &&) = default;

        PredictionBasedQuantizer &operator=(PredictionBasedQuantizer const &) = default;

        PredictionBasedQuantizer &operator=(PredictionBasedQuantizer &&) = default;

        void precompress_data() const {}

        void postcompress_data() const {}

        void predecompress_data() const {}

        void postdecompress_data() const {}

//        void predecompress_block() {}
//
//        void precompress_block() {}


        PredictionBasedQuantizer(T eb, int r) : error_bound(eb),
                                                error_bound_reciprocal(1.0 / eb),
                                                radius(r) {}

        int get_radius() const { return radius; }

        T get_eb() const { return error_bound; }

    };

    template<class T>
    class LinearQuantizer : public PredictionBasedQuantizer<T> {
    public:
        LinearQuantizer(T eb, int r = 32768) : PredictionBasedQuantizer<T>(eb, r) {
        }

        using value_type = T;
        using reference = T &;

        // quantize the data with a prediction value, and returns the quantization index
        int quantize(T data, T pred);

        // quantize the data with a prediction value, and returns the quantization index and the decompressed data
        // int quantize(T data, T pred, T& dec_data);
        int quantize_and_overwrite(T &data, T pred);

        // recover the data using the quantization index
        T recover(T pred, int quant_index);

        void save(unsigned char *&c) const {
            // std::string serialized(sizeof(uint8_t) + sizeof(T) + sizeof(int),0);
            c[0] = 0b00000010;
            c += 1;
            // std::cout << "saving eb = " << this->error_bound << ", unpred_num = "  << unpred.size() << std::endl;
            *reinterpret_cast<T *>(c) = this->error_bound;
            c += sizeof(T);
            *reinterpret_cast<int *>(c) = this->radius;
            c += sizeof(int);
            *reinterpret_cast<size_t *>(c) = unpred.size();
            c += sizeof(size_t);
            memcpy(c, unpred.data(), unpred.size() * sizeof(T));
            c += unpred.size() * sizeof(T);
        };

        void load(const unsigned char *&c, size_t &remaining_length) {
            assert(remaining_length > (sizeof(uint8_t) + sizeof(T) + sizeof(int)));
            c += sizeof(uint8_t);
            remaining_length -= sizeof(uint8_t);
            this->error_bound = *reinterpret_cast<const T *>(c);
            this->error_bound_reciprocal = 1.0 / this->error_bound;
            c += sizeof(T);
            this->radius = *reinterpret_cast<const int *>(c);
            c += sizeof(int);
            size_t unpred_size = *reinterpret_cast<const size_t *>(c);
            c += sizeof(size_t);
            this->unpred = std::vector<T>(reinterpret_cast<const T *>(c), reinterpret_cast<const T *>(c) + unpred_size);
            c += unpred_size * sizeof(T);
            // std::cout << "loading: eb = " << this->error_bound << ", unpred_num = "  << unpred.size() << std::endl;
            // reset index
            index = 0;
        }

        void clear() {
            unpred.clear();
        }

    private:
        std::vector<T> unpred;
        size_t index = 0; // used in decompression only
    };

    template<class T>
    int LinearQuantizer<T>::quantize(T data, T pred) {
        int radius = this->radius;
        // compute quantization index
        int quant_index = (int) ((data - pred) * this->error_bound_reciprocal);
        quant_index = (quant_index > 0) ? (quant_index + 1) / 2 : (quant_index - 1) / 2;
        // shift quantization index, set overbound to 0
        return (quant_index > 0) ? (quant_index < radius ? quant_index + radius : 0) : (quant_index > -radius ? quant_index +
                                                                                                                radius : 0);
    }

//    template<class T>
//    inline int LinearQuantizer<T>::quantize_and_overwrite(T &data, T pred) {
//
//        int quant_index = round((data - pred) * this->error_bound_reciprocal_divided_by_2 + 0.5);
//        int quant = 0;
//        if (abs(quant_index) < this->radius) {
//            T decom = pred + quant_index * this->error_bound_times_2;
//            if (fabs(decom - data) < this->error_bound) {
//                data = decom;
//                quant = quant_index + this->radius;
//            }
//        }
//        if (quant == 0) {
//            unpred.push_back(data);
//        }
//        return quant;
//    }

    template<class T>
    inline int LinearQuantizer<T>::quantize_and_overwrite(T &data, T pred) {

        T diff = data - pred;
        T error_bound = this->error_bound;
        T error_bound_reciprocal = this->error_bound_reciprocal;
        int quant_index_shifted = this->radius + round(diff * error_bound_reciprocal/2.0);
        T decompressed_data = pred + (quant_index_shifted - this->radius)*2 * error_bound;
        if(quant_index_shifted >= 2*this->radius || quant_index_shifted <=0 || fabs(decompressed_data - data) > error_bound) {
            // Handle Unpredictable data
            unpred.push_back(data);
            return 0;// 0 denotes unpredictable data
        }
        data = decompressed_data;
        return quant_index_shifted;
//        int quant_index = (int) (fabs(diff) * this->error_bound_reciprocal) + 1;
//        if (quant_index < this->radius * 2) {
//            quant_index >>= 1;
//            int half_index = quant_index;
//            quant_index <<= 1;
//            int quant_index_shifted;
//            if (diff < 0) {
//                quant_index = -quant_index;
//                quant_index_shifted = this->radius - half_index;
//            } else {
//                quant_index_shifted = this->radius + half_index;
//            }
//            T decompressed_data = pred + quant_index * this->error_bound;
//            if (fabs(decompressed_data - data) > this->error_bound) {
//                unpred.push_back(data);
//                return 0;
//            } else {
//                data = decompressed_data;
//                return quant_index_shifted;
//            }
//        } else {
//            unpred.push_back(data);
//            return 0;
//        }
    }

    template<class T>
    T LinearQuantizer<T>::recover(T pred, int quant_index) {
        if (quant_index) {
            return pred + 2 * (quant_index - this->radius) * this->error_bound;
        } else {
            return unpred[index++];
        }
    }


    // Quantizer for multiple
    template<class T>
    class MultipleErrorBoundsQuantizer : public concepts::QuantizerInterface<T>  {
    public:
        MultipleErrorBoundsQuantizer(std::vector<std::tuple<T,T,T>> eb, int r = 32768):ebs(eb), radius(r){
            init();
        }
        void init(){
            for(int i=0;i<ebs.size();i++){
                T low = std::get<0>(ebs[i]);
                T high = std::get<1>(ebs[i]);
                T eb_cur = std::get<2>(ebs[i]);
                T tmp = (high-low)/(2*eb_cur);
                int quant_num;
                if(tmp> floor(tmp)){
                    quant_num = (int)floor(tmp)+1;
                }else{
                    quant_num = (int)floor(tmp);
                }
//                quant_num = (int)round(tmp);
                std::get<2>(ebs[i]) = (high-low)/(2*quant_num);
                quant_range.push_back(quant_num);
            }
            range_size = ebs.size();
            last_data_range = -1;
            last_data_range = -1;
        }

        void precompress_data() const {}

        void postcompress_data() const {}

        void predecompress_data() const {}

        void postdecompress_data() const {}

        using value_type = T;
        using reference = T &;

        // quantize the data with a prediction value, and returns the quantization index
        int quantize(T data, T pred);
        std::tuple<T, int> quantize_actual(T data, T pred);
        // quantize the data with a prediction value, and returns the quantization index and the decompressed data
        // int quantize(T data, T pred, T& dec_data);
        int quantize_and_overwrite(T &data, T pred);

        // recover the data using the quantization index
        T recover(T pred, int quant_index);

        void save(unsigned char *&c) const {
            // std::string serialized(sizeof(uint8_t) + sizeof(T) + sizeof(int),0);
            c[0] = 0b00000110;
            c += 1;
            // std::cout << "saving eb = " << this->error_bound << ", unpred_num = "  << unpred.size() << std::endl;
            *reinterpret_cast<int *>(c) = ebs.size();
            c += sizeof(int);
            for(int i=0;i<ebs.size();i++){
                *reinterpret_cast<T *>(c) = std::get<0>(ebs[i]);
                c += sizeof(T);
                *reinterpret_cast<T *>(c) = std::get<1>(ebs[i]);
                c += sizeof(T);
                *reinterpret_cast<T *>(c) = std::get<2>(ebs[i]);
                c += sizeof(T);
            }
            *reinterpret_cast<int *>(c) = this->radius;
            c += sizeof(int);
            *reinterpret_cast<size_t *>(c) = unpred.size();
            c += sizeof(size_t);
            memcpy(c, unpred.data(), unpred.size() * sizeof(T));
            c += unpred.size() * sizeof(T);
        };

        void load(const unsigned char *&c, size_t &remaining_length) {
            assert(remaining_length > (sizeof(uint8_t) + sizeof(T) + sizeof(int)));
            c += sizeof(uint8_t);
            remaining_length -= sizeof(uint8_t);
            int size = *reinterpret_cast<const int *>(c);
            c += sizeof(int);
            ebs = std::vector<std::tuple<T,T,T>>(size);
            for(int i=0;i<size;i++){
                T low, high, eb;
                low = *reinterpret_cast<const T *>(c);
                c += sizeof(T);
                high = *reinterpret_cast<const T *>(c);
                c += sizeof(T);
                eb = *reinterpret_cast<const T *>(c);
                c += sizeof(T);
                ebs[i] = std::tuple<T,T,T>(low, high, eb);
            }
//            this->error_bound = *reinterpret_cast<const T *>(c);
//            this->error_bound_reciprocal = 1.0 / this->error_bound;
//            c += sizeof(T);
            this->radius = *reinterpret_cast<const int *>(c);
            c += sizeof(int);
            size_t unpred_size = *reinterpret_cast<const size_t *>(c);
            c += sizeof(size_t);
            this->unpred = std::vector<T>(reinterpret_cast<const T *>(c), reinterpret_cast<const T *>(c) + unpred_size);
            c += unpred_size * sizeof(T);
            // std::cout << "loading: eb = " << this->error_bound << ", unpred_num = "  << unpred.size() << std::endl;
            // reset index
            index = 0;
            init();
        }

        void clear() {
            unpred.clear();
        }
        int get_radius() const { return radius;}

    private:
        std::vector<T> unpred;
        std::vector<int> quant_range;
        size_t index = 0; // used in decompression only
        std::vector<std::tuple<T,T,T>> ebs;
        int radius;
        int last_pred_range, last_data_range;
        int range_size;

        int getErrorBoundIndex(const T& data, bool change_last){
            int cur_range;
            if(data< std::get<1>(ebs[0])){
                cur_range = 0;
            } else if(data>= std::get<0>(ebs[ebs.size()-1])) {
                cur_range = range_size - 1;
            } else {
                if (last_data_range > 0) {
                    cur_range = last_data_range;
                } else {
                    cur_range = 0;
                }
                if (std::get<0>(ebs[cur_range]) <= data) {
                    while (data >= std::get<1>(ebs[cur_range]) && cur_range < range_size - 1) {
                        cur_range++;
                    }
                } else {
                    while (data < std::get<0>(ebs[cur_range]) && cur_range > 0) {
                        cur_range--;
                    }
                }
            }
            if(change_last) {
                last_data_range = cur_range;
            }
            return cur_range;
        }

    };

    template<class T>
    std::tuple<T, int> MultipleErrorBoundsQuantizer<T>::quantize_actual(T data, T pred) {
        // The function returns 0,1 or unshifted quantization value
        T diff = data - pred;
        int data_index = getErrorBoundIndex(data, true);
        int pred_index = getErrorBoundIndex(pred, false);
        int quant_index_shifted;
        T decompressed_data;
        int tmp;
//        assert(pred_index == data_index);
        if(pred_index == data_index){
            T error_bound = std::get<2>(ebs[data_index]);
            T error_bound_reciprocal = 1/error_bound;
            int quant = (int)round(diff *(error_bound_reciprocal*0.5));
            quant_index_shifted = this->radius + quant;
            decompressed_data = pred + quant* 2 * error_bound;
            int tmp_index = getErrorBoundIndex(decompressed_data, false);
            if(pred_index > tmp_index ){
                decompressed_data = std::get<0>(ebs[pred_index]);
            } else if (pred_index < tmp_index) {
                decompressed_data = std::get<1>(ebs[pred_index]);
            }
            auto t1=fabs(decompressed_data - std::get<1>(ebs[tmp_index])),
                t2=std::get<2>(ebs[tmp_index]),
                t3=fabs(decompressed_data - std::get<0>(ebs[tmp_index]));
            if( t1< t2 && fabs(t1-t2)>std::numeric_limits<T>::epsilon()){
                decompressed_data = std::get<1>(ebs[tmp_index]);
            } else if(t3<t2&& fabs(t3-t2)>std::numeric_limits<T>::epsilon()){
                decompressed_data = std::get<0>(ebs[tmp_index]);
            }
        } else if(pred_index > data_index) {
            int quant_value = 0;
            quant_value -= (int)round((pred - std::get<0>(ebs[pred_index]))/(2*std::get<2>(ebs[pred_index])));
            for(int i=pred_index-1; i>=data_index+1;i--){
                quant_value -= quant_range[i];
            }
            tmp = (int)round((std::get<1>(ebs[data_index])-data)/(2*std::get<2>(ebs[data_index])));
            decompressed_data = std::get<1>(ebs[data_index]) - tmp * (2*std::get<2>(ebs[data_index]));
            quant_value -= tmp;
            quant_index_shifted = this->radius + quant_value;
        } else if(pred_index < data_index) {
            int quant_value = 0;
            quant_value += (int)round((std::get<1>(ebs[pred_index]) - pred)/(2*std::get<2>(ebs[pred_index])));
            for(int i=pred_index+1; i<=data_index-1;i++){
                quant_value += quant_range[i];
            }
            tmp = (int)round((data - std::get<0>(ebs[data_index]))/(2*std::get<2>(ebs[data_index])));
            decompressed_data = std::get<0>(ebs[data_index]) + tmp * (2*std::get<2>(ebs[data_index]));
            quant_value += tmp;
            quant_index_shifted = this->radius + quant_value;
        }

        if(quant_index_shifted >= 2*this->radius || quant_index_shifted <=0 || fabs(decompressed_data - data) > std::get<2>(ebs[data_index])) {
            // Handle Unpredictable data
            unpred.push_back(data);
            return std::tuple<T, int>(data, 0); // 0 denotes unpredictable data
        }
        return std::tuple<T, int>(decompressed_data, quant_index_shifted);
    }

    template<class T>
    int MultipleErrorBoundsQuantizer<T>::quantize(T data, T pred){
        auto t = quantize_actual(data, pred);
        return std::get<1>(t);
    }

    template<class T>
    int MultipleErrorBoundsQuantizer<T>::quantize_and_overwrite(T &data, T pred) {
        auto t = quantize_actual(data, pred);
        data = std::get<0>(t);
        return std::get<1>(t);
    }

    template<class T>
    T MultipleErrorBoundsQuantizer<T>::recover(T pred, int quant_index) {
        if(quant_index==0){
            return unpred[index++];
        }
        int pred_index = getErrorBoundIndex(pred, true);
        int actual_quant = quant_index - this->radius;
        int i=pred_index, remaining_quant = actual_quant;
        int tmp;
        T decompressed_data = pred;
        if(actual_quant<0){
            tmp = (int)round((pred - std::get<0>(ebs[pred_index]))/(2*std::get<2>(ebs[pred_index])));
            if(actual_quant + tmp < 0 && pred_index>0){
                remaining_quant += tmp;
                for(i=pred_index-1;i>0;i--){
                    if(remaining_quant+quant_range[i]>=0){
                        decompressed_data = std::get<1>(ebs[i]) + remaining_quant*(2*std::get<2>(ebs[i]));
                        return decompressed_data;
                    }
                    remaining_quant+=quant_range[i];
                }
                // i = 0
                decompressed_data = std::get<1>(ebs[0]) + remaining_quant*(2*std::get<2>(ebs[0]));
                return decompressed_data;
            } else {
                decompressed_data = pred+ remaining_quant*(2*std::get<2>(ebs[pred_index]));
                if(actual_quant+ tmp ==0) {
                    int tmp_index = getErrorBoundIndex(decompressed_data, false);
                    auto t1 = fabs(decompressed_data - std::get<1>(ebs[tmp_index])),
                        t2=std::get<2>(ebs[tmp_index]),
                        t3=fabs(decompressed_data - std::get<0>(ebs[tmp_index]));
                    if ( t1< t2 && fabs(t1-t2)>std::numeric_limits<T>::epsilon()) {
                        decompressed_data = std::get<1>(ebs[tmp_index]);
                    } else if (t3 < t2 && fabs(t3-t2)>std::numeric_limits<T>::epsilon()) {
                        decompressed_data = std::get<0>(ebs[tmp_index]);
                    }
                }
                return decompressed_data;
            }
        } else if(actual_quant == 0){
            int tmp_index = pred_index;
            decompressed_data = pred;
            auto t1 = fabs(decompressed_data - std::get<1>(ebs[tmp_index])),
                    t2=std::get<2>(ebs[tmp_index]),
                    t3=fabs(decompressed_data - std::get<0>(ebs[tmp_index]));
            if ( t1< t2 && fabs(t1-t2)>std::numeric_limits<T>::epsilon()) {
                decompressed_data = std::get<1>(ebs[tmp_index]);
            } else if (t3 < t2 && fabs(t3-t2)>std::numeric_limits<T>::epsilon()) {
                decompressed_data = std::get<0>(ebs[tmp_index]);
            }
            return decompressed_data;
        } else {
            tmp = (int)round((std::get<1>(ebs[pred_index])- pred)/(2*std::get<2>(ebs[pred_index])));
            if(actual_quant - tmp > 0 && pred_index<range_size-1) {
                remaining_quant -= tmp;
                for(i=pred_index+1;i<range_size-1;i++){
                    if(remaining_quant - quant_range[i]<=0){
                        decompressed_data = std::get<0>(ebs[i]) + 2*remaining_quant*std::get<2>(ebs[i]);
                        return decompressed_data;
                    }
                    remaining_quant -=quant_range[i];
                }
                // i = range_size -1
                decompressed_data = std::get<0>(ebs[range_size -1]) + 2*remaining_quant*std::get<2>(ebs[range_size -1]);
                return decompressed_data;
            } else {
                decompressed_data = pred + 2*remaining_quant*std::get<2>(ebs[pred_index]);
                int tmp_index = getErrorBoundIndex(decompressed_data, false);
                auto t1 = fabs(decompressed_data - std::get<1>(ebs[tmp_index])),
                        t2=std::get<2>(ebs[tmp_index]),
                        t3=fabs(decompressed_data - std::get<0>(ebs[tmp_index]));
                if ( t1< t2 && fabs(t1-t2)>std::numeric_limits<T>::epsilon()) {
                    decompressed_data = std::get<1>(ebs[tmp_index]);
                } else if (t3 < t2 && fabs(t3-t2)>std::numeric_limits<T>::epsilon()) {
                    decompressed_data = std::get<0>(ebs[tmp_index]);
                }
                return decompressed_data;
            }
        }
        return decompressed_data;
    }

}
#endif
