#include "Math/math-functions.h"
#include <iostream>
#include <thread>
#include "MyFunctions/my-functions.h"
#include <fmt/core.h>
#include <fmt/color.h>
#include <fmt/ranges.h>
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/async.h"

using namespace sci;

#define _VERIFY_
// #define _CHECK_
#define _PERFORMANCE_

#define MAX_THREADS 4

int party, port = 32000;
int num_threads = 4;
string address = "127.0.0.1";
spdlog::filename_t log_filename;

// 将一维向量转换为列联表，M为列联表的行数，L为列联表的列数
void generate_contingency_table(const std::vector<uint64_t> &hist_vec, int M, int L,
                                std::vector<uint64_t> &R_mat_M, std::vector<uint64_t> &P_mat_L, int bw = 64)
{
    uint64_t mask = my::mask_gen_uint64_t(bw);
    std::vector<uint64_t> R_vec_M(M);
    std::vector<uint64_t> P_vec_L(L);
    for (int j = 0; j < L; j++)
    {
        P_vec_L[j] = 0;
        for (int i = 0; i < M; i++)
        {
            // P_vec_L[j] = (P_vec_L[j] +hist_vec[i * L + j]) & mask;
            P_vec_L[j] = (P_vec_L[j] + hist_vec[i * L + j]);
        }
        P_vec_L[j] = P_vec_L[j] & mask;
    }
    for (int i = 0; i < M; i++)
    {
        R_vec_M[i] = std::reduce(std::execution::par_unseq, hist_vec.begin() + i * L, hist_vec.begin() + (i + 1) * L) & mask;
        std::fill(R_mat_M.begin() + i * L, R_mat_M.begin() + (i + 1) * L, R_vec_M[i]);
        std::copy(P_vec_L.begin(), P_vec_L.end(), P_mat_L.begin() + i * L);
    }
}

// 结果放大了4*scale
ShareA SecureChi2(
    std::vector<ShareA> &sharing_hist_vec,
    uint64_t data_space_left, uint64_t data_space_step,
    int M, int L,
    int bw_hist,
    int scale, int bw_div_res,
    int K = -1)
{
#ifdef _PERFORMANCE_

    std::shared_ptr<spdlog::logger> logger;
    try
    {
        logger = spdlog::basic_logger_st<spdlog::async_factory>("Chi2", log_filename);
    }
    catch (const spdlog::spdlog_ex &ex)
    {
        std::cout << "Log init failed: " << ex.what() << std::endl;
    }
    string party_name = (party == ALICE ? "ALICE" : "BOB");
    auto rounds_begin = iopackArr[0]->get_rounds();

    fmt::print(fg(fmt::color::green), "{}->rounds_begin= {}, bw_hist= {}, scale= {}, bw_div_res= {}\n", party_name, rounds_begin, bw_hist, scale, bw_div_res);
    logger->info("{}->rounds_begin= {}, bw_hist= {}, scale= {}, bw_div_res= {}", party_name, rounds_begin, bw_hist, scale, bw_div_res);

    uint64_t last_rounds = rounds_begin;

    uint64_t total_comm = 0;
    uint64_t thread_comm[num_threads];
    for (int i = 0; i < num_threads; i++)
    {
        thread_comm[i] = iopackArr[i]->get_comm();
    }
    auto start = clock_start();
#endif

    // 获取关键参数

    int len_hist = M * L;
    uint64_t mask_hist = my::mask_gen_uint64_t(bw_hist);
    // uint64_t mask_numeric = my::mask_gen_uint64_t(bw_numeric);
    uint64_t mask_div_res = my::mask_gen_uint64_t(bw_div_res);

    // 初始化计算功能
    auto linearOT_ptr = std::make_unique<LinearOT>(party, iopackArr[0], otpackArr[0]);

    std::vector<ShareA> sharing_R_mat(len_hist);
    std::vector<ShareA> sharing_P_mat(len_hist);

    generate_contingency_table(sharing_hist_vec, M, L, sharing_R_mat, sharing_P_mat);

    std::vector<ShareA> sharing_double_hist_vec = my::append_vectors(sharing_hist_vec, sharing_hist_vec);
    std::vector<ShareA> sharing_RP_vec = my::append_vectors(sharing_R_mat, sharing_P_mat);
    std::vector<ShareA> sharing_uv_vec(len_hist * 2);

    // std::vector<ShareA> sharing_double_hist_scaled_vec(sharing_double_hist_vec);
    // std::transform(std::execution::par_unseq,
    //                sharing_double_hist_scaled_vec.begin(), sharing_double_hist_scaled_vec.end(),
    //                sharing_double_hist_scaled_vec.begin(), [scale](auto e)
    //                { return e << scale; });

    // 多线程除法
    /************** Fork Threads ****************/
    std::thread secure_div[num_threads];
    int chunk_size = len_hist * 2 / num_threads;
    for (int tid = 0; tid < num_threads; ++tid)
    {
        int offset = tid * chunk_size;
        int num_ops;
        if (tid == (num_threads - 1))
        {
            num_ops = len_hist * 2 - offset;
        }
        else
        {
            num_ops = chunk_size;
        }
        // 放大了 scale*2,
        secure_div[tid] = std::thread(
            my::secure_div_thread, party, iopackArr, otpackArr, tid, num_ops, sharing_double_hist_vec.data() + offset,
            sharing_RP_vec.data() + offset, sharing_uv_vec.data() + offset,
            bw_hist, bw_hist, bw_div_res, 0, scale, scale,
            false, true);
    }
    for (int tid = 0; tid < num_threads; ++tid)
    {
        secure_div[tid].join();
    }

#ifdef _PERFORMANCE_
    uint64_t secure_div_rounds = iopackArr[0]->get_rounds() - last_rounds;
    last_rounds = iopackArr[0]->get_rounds();
    fmt::print(fg(fmt::color::white_smoke), "{}->secure_div rounds= {}, num rounds= {}\n", party_name, secure_div_rounds, last_rounds);
    logger->info("{}->secure_div rounds= {}, num rounds= {}", party_name, secure_div_rounds, last_rounds);
#endif

    // 多线程乘法
    std::vector<ShareA> sharing_w_vec(len_hist);
    std::thread secure_hadamard[num_threads];
    chunk_size = len_hist / num_threads;
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
        // 结果放大了 scale*2*2,
        secure_hadamard[tid] = std::thread(
            my::hadamard_product_thread, party, iopackArr, otpackArr, tid,
            num_ops, sharing_uv_vec.data() + offset, sharing_uv_vec.data() + len_hist + offset, sharing_w_vec.data() + offset,
            bw_div_res, bw_div_res, bw_div_res, false, false, MultMode::None, nullptr, nullptr);
    }
    for (int tid = 0; tid < num_threads; ++tid)
    {
        secure_hadamard[tid].join();
    }

#ifdef _PERFORMANCE_
    uint64_t hadamard_product_1_rounds = iopackArr[0]->get_rounds() - last_rounds;
    last_rounds = iopackArr[0]->get_rounds();
    fmt::print(fg(fmt::color::white_smoke), "{}->hadamard_product_1_rounds= {}, num rounds= {}\n", party_name, hadamard_product_1_rounds, last_rounds);
    logger->info("{}->hadamard_product_1_rounds= {}, num rounds= {}", party_name, hadamard_product_1_rounds, last_rounds);
#endif

    //  保持scale*2*2的缩放，求出CHI2. 需要bw_numeric和bw_hist的乘法，开销较大
    ShareA sharing_w = std::reduce(std::execution::par_unseq, sharing_w_vec.begin(), sharing_w_vec.end()) & mask_div_res;
    ShareA sharing_w_minus_1 = (sharing_w - my::shareA_const(party, 1ULL << (scale * 4))) & mask_div_res;

    ShareA sharing_chi2;
    if (K == -1)
    {
        ShareA sharing_K = std::reduce(std::execution::par_unseq, sharing_P_mat.begin(), sharing_P_mat.begin() + L) & mask_hist;
        linearOT_ptr->hadamard_product(1, &sharing_w_minus_1, &sharing_K, &sharing_chi2, bw_div_res, bw_hist, bw_div_res, false, false, MultMode::None, nullptr, nullptr);

#ifdef _PERFORMANCE_
        uint64_t hadamard_product_2_rounds = iopackArr[0]->get_rounds() - last_rounds;
        last_rounds = iopackArr[0]->get_rounds();
        fmt::print(fg(fmt::color::white_smoke), "{}->hadamard_product_2_rounds= {}, num rounds= {}\n", party_name, hadamard_product_2_rounds, last_rounds);
        logger->info("{}->hadamard_product_2_rounds= {}, num rounds= {}", party_name, hadamard_product_2_rounds, last_rounds);
#endif
    }
    else
    {
        sharing_chi2 = sharing_w_minus_1 * K;
    }
    // ShareA sharing_chi2_numeric = sharing_chi2 >> (bw_div_res - bw_numeric);
    // 结果放大了4*scale-(bw_div_res-bw_numeric)

#ifdef _PERFORMANCE_
    /**** Process & Write Benchmarking Data *****/
    long long t = time_from(start);
    for (int i = 0; i < num_threads; i++)
    {
        thread_comm[i] = iopackArr[i]->get_comm() - thread_comm[i];
        total_comm += thread_comm[i];
    }
    fmt::print(fg(fmt::color::white_smoke), "{}->Number of dim/s:\t{}\n", party_name, (double(len_hist) / t) * 1e6);
    fmt::print(fg(fmt::color::white_smoke), "{}->Computing Time\t{} ms\n", party_name, t / (1000.0));
    fmt::print(fg(fmt::color::white_smoke), "{}->Computing Bytes Sent\t{} bytes, {} MB\n", party_name, total_comm, total_comm / (1024.0 * 1024.0));
    logger->info("{}->Number of dim/s:\t{}", party_name, (double(len_hist) / t) * 1e6);
    logger->info("{}->Computing Time\t{} ms", party_name, t / (1000.0));
    logger->info("{}->Computing Bytes Sent\t{} bytes, {} MB", party_name, total_comm, total_comm / (1024.0 * 1024.0));
#endif

    /************** Verification ****************/
#ifdef _VERIFY_

#ifdef _CHECK_
    auto hist_vec = my::reconstruct_uint64_t(party, sci::ALICE, iopackArr, sharing_hist_vec, bw_hist, __LINE__);
    auto R_mat = my::reconstruct_uint64_t(party, sci::ALICE, iopackArr, sharing_R_mat, bw_hist, __LINE__);
    auto P_mat = my::reconstruct_uint64_t(party, sci::ALICE, iopackArr, sharing_P_mat, bw_hist, __LINE__);
    auto double_hist_vec = my::reconstruct_uint64_t(party, sci::ALICE, iopackArr, sharing_double_hist_vec, bw_hist, __LINE__);
    auto RP_vec = my::reconstruct_uint64_t(party, sci::ALICE, iopackArr, sharing_RP_vec, bw_hist, __LINE__);
    auto uv_vec = my::reconstruct_uint64_t(party, sci::ALICE, iopackArr, sharing_uv_vec, bw_div_res, __LINE__);
    auto w_vec = my::reconstruct_uint64_t(party, sci::ALICE, iopackArr, sharing_w_vec, bw_div_res, __LINE__);
#endif

    std::vector<ShareA> sharing_w_minus_1_vec(1, sharing_w_minus_1);
    auto w_minus_1_vec = my::reconstruct_uint64_t(party, sci::ALICE, iopackArr, sharing_w_minus_1_vec, bw_div_res);
    std::vector<ShareA> sharing_chi2_vec(1, sharing_chi2);
    auto chi2 = my::reconstruct_uint64_t(party, sci::ALICE, iopackArr, sharing_chi2_vec, bw_div_res);

    if (party == ALICE)
    {
        fmt::print(fg(fmt::color::blue), "-----------------------Verification in SecChiQuared-----------------------\n");
        fmt::print(fg(fmt::color::blue), "{}->scale= {}\n", party_name, scale);
        fmt::print(fg(fmt::color::blue), "{}->w_minus_1_vec: {}, K (number_of_data) = {}\n", party_name, w_minus_1_vec, K);
        fmt::print(fg(fmt::color::blue), "{}->chi2: {}, scale: {} -> chi2= {}\n", party_name, chi2, 4 * scale, double(chi2[0]) / (1ULL << (4 * scale)));

#ifdef _CHECK_
        fmt::print("{}->hist_vec: {}\n", party_name, hist_vec);
        fmt::print("{}->R_mat: {}\n", party_name, R_mat);
        fmt::print("{}->P_mat: {}\n", party_name, P_mat);
        fmt::print("{}->double_hist_vec: {}\n", party_name, double_hist_vec);
        fmt::print("{}->RP_vec: {}\n", party_name, RP_vec);
        fmt::print("{}->uv_vec: {}\n", party_name, uv_vec);
        fmt::print("{}->w_vec: {}\n", party_name, w_vec);
#endif
    }
#endif

    /******************* Cleanup ****************/
    spdlog::drop_all();
    return sharing_chi2;
}

int main(int argc, char **argv)
{
    int left = -1;
    int step = -1;
    int number_of_data = -1;
    uint64_t generate_min = -1;
    uint64_t generate_max = -1;
    int bw_hist = -1;
    int scale = -1;
    int bw_numeric = -1;
    int bw_div_res = -1;
    int M = -1;
    int L = -1;

    /************* Argument Parsing  ************/
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.arg("po", port, "Port Number");
    amap.arg("nt", num_threads, "Number of threads");
    // amap.arg("lh", len_hist, "Length of histgram vector");
    amap.arg("left", left, "频数向量的左边界");
    amap.arg("step", step, "频数向量的步长");
    amap.arg("sc", scale, "应当是ceil(log2(number_of_data))以上");
    amap.arg("bh", bw_hist, "Bitwidth of histgram");
    amap.arg("bd", bw_div_res, "Bitwidth of res of secdiv, 在chi2中应当是ceil(log2(number_of_data))的4倍以+bw_hist以上");
    amap.arg("min", generate_min, "生成的数据的最小值");
    amap.arg("max", generate_max, "生成的数据的最大值");
    amap.arg("nd", number_of_data, "生成的数据的个数");
    amap.arg("M", M, "M");
    amap.arg("L", L, "L");
    amap.arg("bn", bw_numeric, "Bitwidth of numeric data, 在chi2中应当是ceil(log2(number_of_data))的2倍以上");
    amap.parse(argc, argv);

    my::set_value_int(left, 100);
    my::set_value_int(step, 1);
    my::set_value_int(number_of_data, 500);
    my::set_value_int(M, 5);
    my::set_value_int(L, 20);
    int len_hist = M * L;
    my::set_value_uint64_t(generate_min, left + len_hist * step * 0.02);
    my::set_value_uint64_t(generate_max, left + (len_hist - 1) * step - len_hist * step * 0.02);
    // my::set_value_int(bw_hist, ceil(log2(number_of_data))*2);
    my::set_value_int(bw_hist, ceil(log2(number_of_data))+1);
    int scale_test =  ceil(log2(number_of_data));
    if(bw_hist + 4 * scale_test<61){
        my::set_value_int(scale, ceil(log2(number_of_data)));
    }
    else{
        my::set_value_int(scale, (60-bw_hist)/4);
    }
    my::set_value_int(bw_div_res, bw_hist + 4 * scale);


    my::set_value_int(bw_numeric, ceil(log2(left + (len_hist - 1) * step)));

    log_filename = fmt::format("logs/my_ring_chi2_{}.txt", party);

    assert(num_threads <= MAX_THREADS);

    std::shared_ptr<spdlog::logger> logger;
    try
    {
        logger = spdlog::basic_logger_st<spdlog::async_factory>("main", log_filename);
    }
    catch (const spdlog::spdlog_ex &ex)
    {
        std::cout << "Log init failed: " << ex.what() << std::endl;
    }
    string party_name = (party == ALICE ? "ALICE" : "BOB");
    fmt::print(fg(fmt::color::blue),
               "-------------- {}-> chi2 ----------------------\nnumber_of_data:{}, bw_numeric:{}, len_hist:{}, bw_hist:{}, bw_div_res:{}, scale:{}, left:{}, generate_min:{}, generate_max:{}\n",
               party_name, number_of_data, bw_numeric, len_hist, bw_hist, bw_div_res, scale, left, generate_min, generate_max);
    logger->info("-------------- {}-> chi2 ----------------------\nnumber_of_data:{}, bw_numeric:{}, len_hist:{}, bw_hist:{}, bw_div_res:{}, scale:{}, left:{}, generate_min:{}, generate_max:{}\n",
                 party_name, number_of_data, bw_numeric, len_hist, bw_hist, bw_div_res, scale, left, generate_min, generate_max);

    /********** Setup IO and Base OTs ***********/
    // 每一方有四个线程，每个线程都有一个OTPack
    // for (int i = 0; i < num_threads; i++)
    // {
    //     iopackArr[i] = new IOPack(party, port + i, address);
    //     if (i & 1)
    //     {
    //         otpackArr[i] = new OTPack(iopackArr[i], 3 - party);
    //     }
    //     else
    //     {
    //         otpackArr[i] = new OTPack(iopackArr[i], party);
    //     }
    // }
    my::setup_io_and_ots(party, iopackArr, otpackArr, num_threads, port, address);

    fmt::print("All Base OTs Done\n");

    /************ Generate Test Data ************/

    std::vector<uint64_t> hist_vec_actual_ALICE;
    std::vector<uint64_t> sharing_hist_vec;

    if (party == ALICE)
    {
        auto numeric_data_vec = my::generate_numeric_data_vector(number_of_data, bw_numeric, generate_min, generate_max);
        hist_vec_actual_ALICE = my::numeric_data_to_hist(len_hist, left, step, numeric_data_vec, number_of_data);
        sharing_hist_vec = my::generate_and_send_sharing_vector_plaintext(iopackArr[0], hist_vec_actual_ALICE, bw_hist);
    }
    else if (party == BOB)
    {
        sharing_hist_vec = my::recv_sharing_vector_plaintext_vector(iopackArr[0], len_hist);
    }

    auto chi2 = SecureChi2(sharing_hist_vec, left, step, M, L, bw_hist, scale, bw_div_res, number_of_data);

/************** Verification ****************/
#ifdef _VERIFY_
    if (party == ALICE)
    {
        std::vector<uint64_t> R_mat_M(M * L);
        std::vector<uint64_t> P_mat_L(M * L);
        generate_contingency_table(hist_vec_actual_ALICE, M, L, R_mat_M, P_mat_L);

#ifdef _CHECK_
        fmt::print("-----------------------Verification-----------------------\n");
        // fmt::print("{}->hist_vec_actual_ALICE: {}\n", party_name, hist_vec_actual_ALICE);
        // fmt::print("{}->R_mat_M: {}\n", party_name, R_mat_M);
        // fmt::print("{}->P_mat_L: {}\n", party_name, P_mat_L);

        std::vector<uint64_t> data_accumulate_vec(len_hist);
        data_accumulate_vec[0] = hist_vec_actual_ALICE[0];
        for (int i = 1; i < len_hist; i++)
        {
            data_accumulate_vec[i] = data_accumulate_vec[i - 1] + hist_vec_actual_ALICE[i];
        }
        for (int i = 1; i < len_hist; i++)
        {
            if (i < 25 || i >= len_hist - 10)
            {
                std::cout << "Debug2: data_accumulate_vec" << i << ": " << data_accumulate_vec[i] << std::endl;
            }
        }
#endif
    }
#endif

    for (int i = 0; i < num_threads; i++)
    {
        delete iopackArr[i];
        delete otpackArr[i];
    }
    spdlog::drop_all();
}
