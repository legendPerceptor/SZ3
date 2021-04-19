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
    struct RangeTuple {
        T low, high, eb;
        RangeTuple(T low, T high, T eb):low(low), high(high), eb(eb){}
        RangeTuple():low(0), high(0), eb(0){}
    };

    template<class T, class Q>
    struct tuple2{
        T a;
        Q b;
        tuple2(T a, Q b):a(a), b(b) {}
        tuple2():a(0),b(0){}
    };

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
        MultipleErrorBoundsQuantizer(std::vector<RangeTuple<T>> eb, int r = 32768):ebs(eb), radius(r){
            init();
        }
        void init(){
            for(int i=0;i<ebs.size();i++){
                T low = ebs[i].low;
                T high = ebs[i].high;
                T eb_cur = ebs[i].eb;
                T tmp = (high-low)/(2*eb_cur);
                int quant_num=round(tmp);
                quant_range.push_back(quant_num);
            }
            global_min = ebs[0].low;
            global_max = ebs[ebs.size()-1].high;
            T last_low = ebs[ebs.size()-1].high;
            T last_eb = ebs[ebs.size()-1].eb;
            T first_high = ebs[0].low;
            T first_eb = ebs[0].eb;
            ebs.push_back(RangeTuple<T>(last_low, 10000, last_eb));
            quant_range.push_back((int)(10000-last_low)/(2*last_eb));
            ebs.insert(ebs.begin(),RangeTuple<T>(-10000, first_high, first_eb));
            quant_range.insert(quant_range.begin(), (int)(first_high+10000)/(2*first_eb));
            range_size = ebs.size();
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
        tuple2<T,int> quantize_actual(T data, T pred);
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
                *reinterpret_cast<T *>(c) = ebs[i].low;
                c += sizeof(T);
                *reinterpret_cast<T *>(c) = ebs[i].high;
                c += sizeof(T);
                *reinterpret_cast<T *>(c) = ebs[i].eb;
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
            ebs = std::vector<RangeTuple<T>>(size);
            for(int i=0;i<size;i++){
                T low, high, eb;
                low = *reinterpret_cast<const T *>(c);
                c += sizeof(T);
                high = *reinterpret_cast<const T *>(c);
                c += sizeof(T);
                eb = *reinterpret_cast<const T *>(c);
                c += sizeof(T);
                ebs[i] = RangeTuple(low, high, eb);
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
        std::vector<RangeTuple<T>> ebs;
        T global_min, global_max;
        int radius;
        int last_data_range;
        int range_size;
        const T EPSILON = 0.001;

        int getErrorBoundIndex(const T& data, bool change_last){
            int cur_range;
            if(data< ebs[0].high){
                cur_range = 0;
            } else if(data>= ebs[ebs.size()-1].low) {
                cur_range = range_size - 1;
            } else {
                if (last_data_range > 0) {
                    cur_range = last_data_range;
                } else {
                    cur_range = 0;
                }
                if (ebs[cur_range].low <= data) {
                    while (data >= ebs[cur_range].high && cur_range < range_size - 1) {
                        cur_range++;
                    }
                } else {
                    while (data < ebs[cur_range].low && cur_range > 0) {
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
    tuple2<T, int> MultipleErrorBoundsQuantizer<T>::quantize_actual(T data, T pred) {
        // The function returns 0,1 or unshifted quantization value
        int tp = round(pred/EPSILON);
        pred = (T)tp * EPSILON;
        int pred_index = getErrorBoundIndex(pred, true);
        if(fabs(ebs[pred_index].high-pred)<EPSILON){
            pred = ebs[pred_index].high;
            pred_index +=1;
        } else if(fabs(ebs[pred_index].low-pred)<EPSILON) {
            pred = ebs[pred_index].low;
        }
        T diff = data - pred;
        int data_index = getErrorBoundIndex(data, true);
//        int pred_index = getErrorBoundIndex(pred, false);
        int quant_index_shifted;
        T decompressed_data;
        int tmp;
//        assert(pred_index == data_index);
        if(pred_index == data_index){
            T error_bound = ebs[data_index].eb;
            T error_bound_reciprocal = 1/error_bound;
            int quant = (int)round(diff *(error_bound_reciprocal*0.5));

            decompressed_data = pred + quant* 2 * error_bound;
//            int tmp_index = getErrorBoundIndex(decompressed_data, false);
            if(fabs(ebs[pred_index].low-decompressed_data)<EPSILON){
                decompressed_data = ebs[pred_index].low;
            }else if(fabs(ebs[pred_index].high-decompressed_data)<EPSILON){
                decompressed_data = ebs[pred_index].high;
            }else if(ebs[pred_index].low > decompressed_data){
                decompressed_data = ebs[pred_index].low + ebs[pred_index].eb;
            } else if (ebs[pred_index].high < decompressed_data) {
                decompressed_data = ebs[pred_index].high - ebs[pred_index].eb;
            }
            quant_index_shifted = this->radius + quant;
        } else if(pred_index > data_index) {
            int quant_value = 0;
            T t_diff = pred - ebs[pred_index].low;
            int t1 = (int)((t_diff)/ebs[pred_index].eb);
            int t2 = (int)round((t_diff)/ebs[pred_index].eb);
            if(t1==t2) {
                quant_value -= t1%2==0? t1/2: (t1+1)/2;
            }else{
                if(fabs((T)t2*ebs[pred_index].eb-t_diff)<EPSILON){
                    quant_value -= t2%2==0 ? t2/2 : (t2+1)/2;
                }else{
                    quant_value -= t1%2==0 ? t1/2 : (t1+1)/2;
                }
            }
            for(int i=pred_index-1; i>=data_index+1;i--){
                quant_value -= quant_range[i];
            }
            t_diff = ebs[data_index].high-ebs[data_index].eb-data;
            t1 = (int) ((t_diff/ebs[data_index].eb));
            t2 = (int) round(t_diff/ebs[data_index].eb);
            if(t1==t2){
                tmp= t1%2==0? t1/2 : (t1+1)/2;
            }else{
                if(fabs((T)t2*ebs[data_index].eb-t_diff)<EPSILON){
                    tmp = t2%2==0 ? t2/2 : (t2+1)/2;
                }else{
                    tmp = t1%2==0 ? t1/2 : (t1+1)/2;
                }
            }
//            tmp = (int)(round((ebs[data_index].high-ebs[data_index].eb-data)/ebs[data_index].eb-0.49999+std::numeric_limits<T>::epsilon())+1)/2;
            decompressed_data = ebs[data_index].high-ebs[data_index].eb - tmp * (2*ebs[data_index].eb);
            if(fabs(decompressed_data-ebs[data_index].low)<EPSILON){
                decompressed_data = ebs[data_index].low;
            }else if(decompressed_data < ebs[data_index].low ){
                decompressed_data = ebs[data_index].low + ebs[data_index].eb;
                if(tmp==quant_range[data_index]){
                    tmp-=1;
                }
            }
            quant_value -= tmp+1;//Add one more quantization value
            quant_index_shifted = this->radius + quant_value;
        } else if(pred_index < data_index) {
            int quant_value = 0;
            T t_diff = ebs[pred_index].high - pred;
            int t1 = (int)((t_diff)/ebs[pred_index].eb);
            int t2 = (int)round((t_diff)/ebs[pred_index].eb);
            if(t1==t2) {
                quant_value += t1%2==0? t1/2: (t1+1)/2;
            }else{
                if(fabs((T)t2*ebs[pred_index].eb-t_diff)<EPSILON){
                    quant_value += t2%2==0 ? t2/2 : (t2+1)/2;
                }else{
                    quant_value += t1%2==0 ? t1/2 : (t1+1)/2;
                }
            }
//            quant_value += (int)(round((ebs[pred_index].high - pred)/ebs[pred_index].eb-0.49999+std::numeric_limits<T>::epsilon())+1)/2;
            for(int i=pred_index+1; i<=data_index-1;i++){
                quant_value += quant_range[i];
            }
            t_diff = data - (ebs[data_index].low+ebs[data_index].eb);
            t1 = (int) (t_diff/ebs[data_index].eb);
            t2 = (int) round(t_diff/ebs[data_index].eb);
            if(t1==t2){
                tmp= t1%2==0? t1/2 : (t1+1)/2;
            }else{
                if(fabs((T)t2*ebs[data_index].eb-t_diff)<EPSILON){
                    tmp = t2%2==0 ? t2/2 : (t2+1)/2;
                }else{
                    tmp = t1%2==0 ? t1/2 : (t1+1)/2;
                }
            }
//            tmp = (int)(round((data - (ebs[data_index].low+ebs[data_index].eb))/ebs[data_index].eb-0.49999+std::numeric_limits<T>::epsilon())+1)/2;
            decompressed_data = ebs[data_index].low+ebs[data_index].eb + tmp * (2*ebs[data_index].eb);
            if(fabs(decompressed_data-ebs[data_index].high)<EPSILON){
                decompressed_data = ebs[data_index].high;
            }else if(decompressed_data > ebs[data_index].high){
                decompressed_data = ebs[data_index].high - ebs[data_index].eb;
                if(tmp==quant_range[data_index]) {
                    tmp-=1;
                }
            }
            quant_value += tmp+1;
            quant_index_shifted = this->radius + quant_value;
        }

        if(quant_index_shifted >= 2*this->radius || quant_index_shifted <=0 || fabs(decompressed_data - data) > ebs[data_index].eb) {
            // Handle Unpredictable data
            unpred.push_back(data);
            return tuple2<T, int>(data, 0); // 0 denotes unpredictable data
        }
        return tuple2<T, int>(decompressed_data, quant_index_shifted);
    }

    template<class T>
    int MultipleErrorBoundsQuantizer<T>::quantize(T data, T pred){
        data = fmin(fmax(global_min, data), global_max);
//        pred = fmin(fmax(global_min, pred), global_max);
        auto t = quantize_actual(data, pred);
        return t.b;
    }

    template<class T>
    int MultipleErrorBoundsQuantizer<T>::quantize_and_overwrite(T &data, T pred) {
        data = fmin(fmax(global_min, data), global_max);
//        pred = fmin(fmax(global_min, pred), global_max);
        auto t = quantize_actual(data, pred);
        data = t.a;
        return t.b;
    }

    template<class T>
    T MultipleErrorBoundsQuantizer<T>::recover(T pred, int quant_index) {
//        pred = fmin(fmax(global_min, pred), global_max);
        if(quant_index==0){
            return unpred[index++];
        }
        int tp = round(pred/EPSILON);
        pred = (T)tp * EPSILON;
        int pred_index = getErrorBoundIndex(pred, true);
        if(fabs(ebs[pred_index].high-pred)<EPSILON){
            pred = ebs[pred_index].high;
            pred_index +=1;
        } else if(fabs(ebs[pred_index].low-pred)<EPSILON) {
            pred = ebs[pred_index].low;
        }
        int actual_quant = quant_index - this->radius;
        int i, remaining_quant = actual_quant;
        int tmp;
        T decompressed_data = pred;
        if(actual_quant<0){
            T t_diff = pred - ebs[pred_index].low;
            int t1 = (int) (t_diff/ebs[pred_index].eb);
            int t2 = (int) round(t_diff/ebs[pred_index].eb);
            if(t1==t2){
                tmp= t1%2==0? t1/2 : (t1+1)/2;
            }else{
                if(fabs((T)t2*ebs[pred_index].eb-t_diff)<EPSILON){
                    tmp = t2%2==0 ? t2/2 : (t2+1)/2;
                }else{
                    tmp = t1%2==0 ? t1/2 : (t1+1)/2;
                }
            }
//            tmp = (int)(round((pred - ebs[pred_index].low)/ebs[pred_index].eb-0.49999+std::numeric_limits<T>::epsilon())+1)/2;
            if(actual_quant + tmp < 0 && pred_index>0){
                remaining_quant += tmp;
                for(i=pred_index-1;i>0;i--){
                    if(remaining_quant+quant_range[i]>=0){
                        decompressed_data = ebs[i].high-ebs[i].eb + (remaining_quant+1)*(2*ebs[i].eb);
                        if(decompressed_data<ebs[i].low){
                            decompressed_data = ebs[i].low+ebs[i].eb;
                        }
                        return decompressed_data;
                    }
                    remaining_quant+=quant_range[i];
                }
                // i = 0
                decompressed_data = ebs[0].high-ebs[0].eb + (remaining_quant+1)*(2*ebs[0].eb);
                return decompressed_data;
            } else { //data and pred are in the same range
                decompressed_data = pred+ remaining_quant*(2*ebs[pred_index].eb);
                if(fabs(decompressed_data-ebs[pred_index].low)< EPSILON){
                    decompressed_data = ebs[pred_index].low;
                }else if(decompressed_data < ebs[pred_index].low){
                    decompressed_data = ebs[pred_index].low + ebs[pred_index].eb;
                }
//                if(actual_quant+ tmp ==0) {
//                    int tmp_index = getErrorBoundIndex(pred, false);
//                    if(decompressed_data < ebs[tmp_index].low) {
//                        decompressed_data = ebs[tmp_index].low;
//                    }
//                }
                return decompressed_data;
            }
        } else if(actual_quant == 0){
//            int tmp_index = pred_index;
            decompressed_data = pred;
//            auto t1 = fabs(decompressed_data - ebs[tmp_index].high),
//                    t2=ebs[tmp_index].eb,
//                    t3=fabs(decompressed_data -ebs[tmp_index].low);
//            if ( t1< t2 && fabs(t1-t2)>std::numeric_limits<T>::epsilon()) {
//                decompressed_data = ebs[tmp_index].high;
//            } else if (t3 < t2 && fabs(t3-t2)>std::numeric_limits<T>::epsilon()) {
//                decompressed_data = ebs[tmp_index].low;
//            }
            return decompressed_data;
        } else {
            T t_diff = ebs[pred_index].high- pred;
            int t1 = (int) (t_diff/ebs[pred_index].eb);
            int t2 = (int) round(t_diff/ebs[pred_index].eb);
            if(t1==t2){
                tmp= t1%2==0? t1/2 : (t1+1)/2;
            }else{
                if(fabs((T)t2*ebs[pred_index].eb-t_diff)<EPSILON){
                    tmp = t2%2==0 ? t2/2 : (t2+1)/2;
                }else{
                    tmp = t1%2==0 ? t1/2 : (t1+1)/2;
                }
            }
//            tmp = (int)(round((ebs[pred_index].high- pred)/ebs[pred_index].eb-0.49999+std::numeric_limits<T>::epsilon())+1)/2;
            if(actual_quant - tmp > 0 && pred_index<range_size-1) {
                remaining_quant -= tmp;
                for(i=pred_index+1;i<range_size-1;i++){
                    if(remaining_quant - quant_range[i]<=0){
                        decompressed_data = ebs[i].low + ebs[i].eb + 2*(remaining_quant-1)*ebs[i].eb;
                        if(decompressed_data>ebs[i].high){
                            decompressed_data = ebs[i].high - ebs[i].eb;
                        }
                        return decompressed_data;
                    }
                    remaining_quant -=quant_range[i];
                }
                // i = range_size -1
                decompressed_data = ebs[range_size -1].low+ebs[range_size-1].eb + 2*(remaining_quant-1)*ebs[range_size -1].eb;
                return decompressed_data;
            } else { // data and pred are in the same range
                decompressed_data = pred + 2*remaining_quant*ebs[pred_index].eb;
                if(fabs(decompressed_data-ebs[pred_index].high)< EPSILON){
                    decompressed_data = ebs[pred_index].high;
                }else if(decompressed_data > ebs[pred_index].high){
                    decompressed_data = ebs[pred_index].high-ebs[pred_index].eb;
                }
//                if(actual_quant- tmp ==0) {
//                    int tmp_index = getErrorBoundIndex(pred, false);
//                    if(decompressed_data>ebs[tmp_index].high) {
//                        decompressed_data = ebs[tmp_index].high-ebs[tmp_index].eb;
//                    }
//                }
                return decompressed_data;
            }
        }
        return decompressed_data;
    }

}
#endif
