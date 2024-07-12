#include "Math/math-functions.h"
#include <iostream>
#include <thread>
#include "MyFunctions/my-functions.h"
#include <vector>
#include <algorithm>

#include <fmt/core.h>
#include <fmt/color.h>

using namespace sci;

// #define _VERIFY_

#define MAX_THREADS 4

int party, port = 32000;
int num_threads = 4;
string address = "127.0.0.1";
// sci::IOPack *iopackArr[MAX_THREADS];
// sci::OTPack *otpackArr[MAX_THREADS];

/*
Q(p, vector x) = (1-gamma) * x_index + gamma * x_{index+1}
index = floor((n-1) * p) + 1
gamma = (n-1) * p - floor((n-1) * p)
n is the length of vector x
x begins from 1
*/
struct QuantileInfo
{
    double percentile;
    int num_data;
    int index;
    double gamma;
};

struct QuantileArithmetic
{
    uint64_t left;
    uint64_t right;
};

QuantileArithmetic quantile_actual;

QuantileInfo quantile_number(double percentile, int num_data)
{
    QuantileInfo qinfo{};
    qinfo.percentile = percentile;
    qinfo.num_data = num_data;
    qinfo.index = floor((num_data - 1) * percentile) + 1;
    qinfo.gamma = (num_data - 1) * percentile - floor((num_data - 1) * percentile);
    return qinfo;
}

QuantileArithmetic hist_quantile(
    double percentile, int num_data,
    std::vector<ShareA> &sharing_hist_vec,
    int bw_hist,
    int len_hist,
    uint64_t data_space_left,
    uint64_t data_space_step,
    int bw_numeric)
{
    QuantileInfo qinfo = quantile_number(percentile, num_data);
    uint64_t compare_value = qinfo.index;

    int bw_compare_result = bw_numeric;

#ifdef _VERIFY_
    std::cout << "Debug4: bitlength_compare_result=" << bw_compare_result << std::endl;
#endif

    uint64_t total_comm = 0;
    uint64_t thread_comm[num_threads];
    for (int i = 0; i < num_threads; i++)
    {
        thread_comm[i] = iopackArr[i]->get_comm();
    }
    auto start = clock_start();

    uint64_t mask_compare_result = my::mask_gen_uint64_t(bw_compare_result);
    uint64_t mask_numeric = my::mask_gen_uint64_t(bw_numeric);
    uint64_t mask_hist = my::mask_gen_uint64_t(bw_hist);

    // TODO 修改1
    auto sharing_v_vec = std::vector<ShareA>(len_hist);
    sharing_v_vec[0] = sharing_hist_vec[0];
    for (int i = 1; i < len_hist; i++)
    {
        sharing_v_vec[i] = (sharing_v_vec[i - 1] + sharing_hist_vec[i]) & mask_hist;
    }
    auto less_vec_B = std::vector<ShareB>(len_hist); // 比较大小返回值向量
    auto less_vec_A = std::vector<ShareA>(len_hist);
    auto eq_vec_B = std::vector<ShareB>(len_hist); // 比较大小返回值向量
    auto eq_vec_A = std::vector<ShareA>(len_hist);

    /************** Fork Threads ****************/
    std::thread compare_arithmetic_threads[num_threads];
    int chunk_size = len_hist / num_threads;
    for (int tid = 0; tid < num_threads; ++tid)
    {
        int offset = tid * chunk_size;
        int num_ops;
        if (tid == (num_threads - 1))
        {
            num_ops = len_hist - offset;
        }
        else
        {
            num_ops = chunk_size;
        }
        compare_arithmetic_threads[tid] = std::thread(
            my::compare_arithmetic_thread, LESSwithEQUAL, party, iopackArr, otpackArr, tid,
            sharing_v_vec.data() + offset, compare_value, num_ops, bw_hist, bw_compare_result,
            less_vec_B.data() + offset, eq_vec_B.data() + offset, less_vec_A.data() + offset, eq_vec_A.data() + offset);
    }
    for (int i = 0; i < num_threads; ++i)
    {
        compare_arithmetic_threads[i].join();
    }

    uint64_t count_less_eq_A_res[2] = {0}; // 0是less，1是equal

    count_less_eq_A_res[0] = std::reduce(std::execution::par_unseq,less_vec_A.begin(), less_vec_A.begin() + len_hist, 0) & mask_compare_result;
    count_less_eq_A_res[1] = std::reduce(std::execution::par_unseq,eq_vec_A.begin(), eq_vec_A.begin() + len_hist, 0) & mask_compare_result;

    // count_less_eq_A_res[0] = std::accumulate(less_vec_A.begin(), less_vec_A.begin() + len_hist, 0) & mask_compare_result;
    // count_less_eq_A_res[1] = std::accumulate(eq_vec_A.begin(), eq_vec_A.begin() + len_hist, 0) & mask_compare_result;


    // long long t = time_from(start);
    // for (int i = 0; i < num_threads; i++)
    // {
    //     thread_comm[i] = iopackArr[i]->get_comm() - thread_comm[i];
    //     total_comm += thread_comm[i];
    // }

    uint64_t count_less_eq_A_num[2]; // 1的数量是max结果的index的sharing

    std::unique_ptr<XTProtocol> xt_ptr = std::make_unique<XTProtocol>(party, iopackArr[0], otpackArr[0]);
    xt_ptr->z_extend(2, count_less_eq_A_res, count_less_eq_A_num, bw_compare_result, bw_numeric);

    // std::cout<<party<<"  >>>"<<count_equal<<"   >>>"<<count_equal_arithmetic<<std::endl;
    // std::cout<<mask_compare_result<<"  >>>"<<mask_numeric<<std::endl;

    // TODO 修改2
    QuantileArithmetic quantile_arithmetic{};

    if (party == ALICE)
    {
        quantile_arithmetic.left = (data_space_left + count_less_eq_A_num[0] * data_space_step) & mask_numeric;
        quantile_arithmetic.right = (data_space_left + (count_less_eq_A_num[0] + count_less_eq_A_num[1]) * data_space_step) & mask_numeric;
    }
    else
    {
        quantile_arithmetic.left = (count_less_eq_A_num[0] * data_space_step) & mask_numeric;
        quantile_arithmetic.right = ((count_less_eq_A_num[0] + count_less_eq_A_num[1]) * data_space_step) & mask_numeric;
    }

    long long t = time_from(start);
    for (int i = 0; i < num_threads; i++)
    {
        thread_comm[i] = iopackArr[i]->get_comm() - thread_comm[i];
        total_comm += thread_comm[i];
    }

    /************** Verification ****************/
#ifdef _VERIFY_
    {
        if (party == ALICE)
        {
            iopackArr[0]->io->send_data(less_vec_B.data(), len_hist * sizeof(uint8_t));
            iopackArr[0]->io->send_data(less_vec_A.data(), len_hist * sizeof(uint64_t));
            iopackArr[0]->io->send_data(eq_vec_B.data(), len_hist * sizeof(uint8_t));
            iopackArr[0]->io->send_data(eq_vec_A.data(), len_hist * sizeof(uint64_t));
            iopackArr[0]->io->send_data(sharing_v_vec.data(), len_hist * sizeof(uint64_t));
        }
        else if (party == BOB)
        {
            std::cout << "------------------------Verifying...------------------------" << std::endl;

            auto res_cmp_ptr0 = std::make_unique<uint8_t[]>(len_hist);
            auto res_cmp_arithmetic_ptr0 = std::make_unique<uint64_t[]>(len_hist);
            auto res_eq_ptr0 = std::make_unique<uint8_t[]>(len_hist);
            auto res_eq_arithmetic_ptr0 = std::make_unique<uint64_t[]>(len_hist);
            auto frequency_V_ptr0 = std::make_unique<uint64_t[]>(len_hist);

            iopackArr[0]->io->recv_data(res_cmp_ptr0.get(), len_hist * sizeof(uint8_t));
            iopackArr[0]->io->recv_data(res_cmp_arithmetic_ptr0.get(), len_hist * sizeof(uint64_t));
            iopackArr[0]->io->recv_data(res_eq_ptr0.get(), len_hist * sizeof(uint8_t));
            iopackArr[0]->io->recv_data(res_eq_arithmetic_ptr0.get(), len_hist * sizeof(uint64_t));
            iopackArr[0]->io->recv_data(frequency_V_ptr0.get(), len_hist * sizeof(uint64_t));

            for (int i = 0; i < len_hist; i++)
            {
                uint8_t Res_cmp = unsigned_val(less_vec_B[i] + res_cmp_ptr0[i], 1);
                uint8_t Res_eq = unsigned_val(eq_vec_B[i] + res_eq_ptr0[i], 1);
                uint64_t Data = unsigned_val(sharing_v_vec[i] + frequency_V_ptr0[i], bw_hist);
                if ((i < 150) || (i >= len_hist - 10))
                {
                    std::cout << "Debug5: " << i << ": "
                              << "Res_cmp: " << int(Res_cmp) << ", " << int(res_cmp_ptr0[i]) << ", " << int(less_vec_B[i]) << "\t"
                              << "Res_cmp_arithmetic: " << res_cmp_arithmetic_ptr0[i] << "  " << less_vec_A[i] << "\t"
                              << "Res_eq: " << int(Res_eq) << ", " << int(res_eq_ptr0[i]) << ", " << int(eq_vec_B[i]) << "\t"
                              << "Res_eq_arithmetic: " << res_eq_arithmetic_ptr0[i] << "  " << eq_vec_A[i] << "\t"
                              << "frequency_V: " << Data << std::endl;
                }

                uint8_t expected_Res_cmp = 0;
                // TODO 修改3
                // todo 需要修改计算方法
                if (i < (quantile_actual.left - data_space_left) / data_space_step)
                {
                    expected_Res_cmp = 1;
                }
                if (Res_cmp != expected_Res_cmp)
                {
                    std::cout << "assert(Res_cmp == expected_Res_cmp);  " << int(i) << " " << int(Res_cmp) << std::endl;
                }
                assert(Res_cmp == expected_Res_cmp);
            }
            std::cout << "------------------------Verification passed!------------------------" << std::endl;
        }
    }
#endif

    /**** Process & Write Benchmarking Data *****/
    string party_name = (party == ALICE ? "ALICE" : "BOB");
    fmt::print(fg(fmt::color::white_smoke), "{}->Number of dim/s:\t{}\n", party_name, (double(len_hist) / t) * 1e6);
    fmt::print(fg(fmt::color::white_smoke), "{}->Computing Time\t{} ms\n", party_name, t / (1000.0));
    fmt::print(fg(fmt::color::white_smoke), "{}->Computing Bytes Sent\t{} bytes\n", party_name, total_comm);

    /******************* Cleanup ****************/
    return quantile_arithmetic;
}

int main(int argc, char **argv)
{
    int len_hist = -1;
    int left = -1;
    int step = -1;
    int number_of_data = -1;
    uint64_t generate_min = -1;
    uint64_t generate_max = -1;
    int bw_hist = -1;
    int bw_numeric = -1;
    double percentile = -1.0;

    /************* Argument Parsing  ************/
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.arg("po", port, "Port Number");
    amap.arg("nt", num_threads, "Number of threads");
    amap.arg("lh", len_hist, "Length of histgram vector");
    amap.arg("bh", bw_hist, "Bitwidth of histgram");
    amap.arg("left", left, "频数向量的左边界");
    amap.arg("step", step, "频数向量的步长");
    amap.arg("bn", bw_numeric, "Bitwidth of numeric data");
    amap.arg("min", generate_min, "生成的数据的最小值");
    amap.arg("max", generate_max, "生成的数据的最大值");
    amap.arg("nd", number_of_data, "生成的数据的个数");
    amap.arg("pe", percentile, "percentile");
    amap.parse(argc, argv);

    my::set_value_int(len_hist, 1000);
    my::set_value_int(left, 100);
    my::set_value_int(step, 1);
    my::set_value_int(number_of_data, 500);
    my::set_value_uint64_t(generate_min, left + left * 0.02);
    my::set_value_uint64_t(generate_max, left + (len_hist - 1) * step - step * 0.02);
    my::set_value_int(bw_hist, 16);
    my::set_value_int(bw_numeric, 14);
    my::set_value_double(percentile, 0.1);

    assert(num_threads <= MAX_THREADS);

    fmt::print(fg(fmt::color::blue), "len_hist:{}, left:{}, generate_min:{}, generate_max:{}, percentile:{}\n", len_hist, left, generate_min, generate_max, percentile);

    /********** Setup IO and Base OTs ***********/
    // 每一方有四个线程，每个线程都有一个OTPack
    for (int i = 0; i < num_threads; i++)
    {
        iopackArr[i] = new IOPack(party, port + i, address);
        if (i & 1)
        {
            otpackArr[i] = new OTPack(iopackArr[i], 3 - party);
        }
        else
        {
            otpackArr[i] = new OTPack(iopackArr[i], party);
        }
    }
    fmt::print("All Base OTs Done\n");

    /************ Generate Test Data ************/
    std::vector<uint64_t> sharing_hist_vec;
    // TODO修改4

    QuantileInfo qinfo = quantile_number(percentile, number_of_data);
    if (party == ALICE)
    {
#ifdef _VERIFY_
        printf("Debug3: qinfo{ percentile:%f, num_data:%d, index:%d, gamma:%f}\n", qinfo.percentile, qinfo.num_data, qinfo.index, qinfo.gamma);
#endif

        auto numeric_data_vec = my::generate_numeric_data_vector(number_of_data, bw_numeric, generate_min, generate_max);
        auto hist_vec = my::numeric_data_to_hist(len_hist, left, step, numeric_data_vec, number_of_data);
        std::vector<uint64_t> data_accumulate_vec(len_hist);
        data_accumulate_vec[0] = hist_vec[0];

        //--------------------------------------------------------
        for (int i = 1; i < len_hist; i++)
        {
            data_accumulate_vec[i] = data_accumulate_vec[i - 1] + hist_vec[i];
#ifdef _VERIFY_
            if (i < 25 || i >= len_hist - 10)
            {
                std::cout << "Debug2: data_accumulate_vec" << i << ": " << data_accumulate_vec[i] << std::endl;
            }
#endif
        }

        // // 原先明文找index的方法
        // int count = 0;
        // for (int i = 0; i < len_hist; i++)
        // {
        //     if (data_accumulate_ptr[i] < qinfo.index)
        //     {
        //         count++;
        //     }
        //     else
        //         break;
        // }
        // quantile_actual_left = count * step + left;
        // std::cout<<"Debug1: count:" << count << " quantile_actual_left: " << quantile_actual_left<<" index " << data_accumulate_ptr[count] << std::endl;

        vector<uint64_t> data_vec(numeric_data_vec);
        std::sort(data_vec.begin(), data_vec.end());
        quantile_actual.left = data_vec[qinfo.index - 1];
        quantile_actual.right = data_vec[qinfo.index];
        fmt::print(fg(fmt::color::blue), "Line:{}, quantile_actual_left={}, quantile_actual_right={}\n", __LINE__, quantile_actual.left, quantile_actual.right);
        iopackArr[0]->io->send_data(&quantile_actual, sizeof(QuantileArithmetic));
        //--------------------------------------------------------

        sharing_hist_vec = my::generate_and_send_sharing_vector_plaintext(iopackArr[0], hist_vec, bw_hist);
    }
    else if (party == BOB)
    {
        //--------------------------------------------------------
        iopackArr[0]->io->recv_data(&quantile_actual, sizeof(QuantileArithmetic));
        //--------------------------------------------------------

        sharing_hist_vec = my::recv_sharing_vector_plaintext_vector(iopackArr[0], len_hist);
    }

    // uint64_t count_equal = hist_quantile(p, number_of_data, sharing_hist_ptr, bw_frequency, len_hist, left, step, bw_arithmetic);
    QuantileArithmetic quantile_arithmetic = hist_quantile(percentile, number_of_data, sharing_hist_vec, bw_hist, len_hist, left, step, bw_numeric);

    /************** Verification ****************/
    if (party == ALICE)
    {
        iopackArr[0]->io->send_data(&quantile_arithmetic, sizeof(QuantileArithmetic));
    }
    else if (party == BOB)
    {
        QuantileArithmetic quantile_arithmetic0;
        iopackArr[0]->io->recv_data(&quantile_arithmetic0, sizeof(QuantileArithmetic));

        // TODO 修改5
        uint64_t quantile_arithmetic_left = unsigned_val(quantile_arithmetic.left + quantile_arithmetic0.left,
                                                         bw_numeric);
        uint64_t quantile_arithmetic_right = unsigned_val(quantile_arithmetic.right + quantile_arithmetic0.right,
                                                          bw_numeric);
        // #ifdef _VERIFY_
        fmt::print(fg(fmt::color::crimson), "Line:{}, Alice:{} Bob:{} quantile_arithmetic_left={}\n", __LINE__, quantile_arithmetic0.left, quantile_arithmetic.left, quantile_arithmetic_left);
        fmt::print(fg(fmt::color::crimson), "Line:{}, Alice:{} Bob:{} quantile_arithmetic_right={}\n", __LINE__, quantile_arithmetic0.right, quantile_arithmetic.right, quantile_arithmetic_right);

        // #endif
        assert((quantile_actual.left - left) / step == (quantile_arithmetic_left - left) / step);
        assert((quantile_actual.right - left) / step == (quantile_arithmetic_right - left) / step);
    }

    /******************* Cleanup ****************/
    for (int i = 0; i < num_threads; i++)
    {
        delete iopackArr[i];
        delete otpackArr[i];
    }
}
