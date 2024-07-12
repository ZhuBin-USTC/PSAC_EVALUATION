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

// 结果放大了2, bw_hist=bw_numeric时，最后乘法不需要扩展位宽，交互更少。
std::vector<ShareA> Secure_U_Test(
    std::vector<ShareA> &sharing_hist_X_vec,
    std::vector<ShareA> &sharing_hist_Y_vec,
    size_t len_hist,
    int bw_hist)
{
#ifdef _PERFORMANCE_
    std::shared_ptr<spdlog::logger> logger;
    try
    {
        logger = spdlog::basic_logger_st<spdlog::async_factory>("U", log_filename);
    }
    catch (const spdlog::spdlog_ex &ex)
    {
        std::cout << "Log init failed: " << ex.what() << std::endl;
    }
    string party_name = (party == ALICE ? "ALICE" : "BOB");
    auto rounds_begin = iopackArr[0]->get_rounds();
    fmt::print(fg(fmt::color::green), "{}->rounds_begin= {}, bw_hist= {}\n", party_name, rounds_begin, bw_hist);
    logger->info("{}->rounds_begin= {}, bw_hist= {}", party_name, rounds_begin, bw_hist);
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
    uint64_t mask_hist = my::mask_gen_uint64_t(bw_hist);

    std::vector<ShareA> sharing_v_X_vec(len_hist);
    std::vector<ShareA> sharing_v_Y_vec(len_hist);
    sharing_v_X_vec[len_hist - 1] = 0;
    sharing_v_Y_vec[len_hist - 1] = 0;
    for (int i = len_hist - 2; i >= 0; i--)
    {
        sharing_v_X_vec[i] = (sharing_hist_X_vec[i + 1] + sharing_v_X_vec[i + 1]) & mask_hist;
        sharing_v_Y_vec[i] = (sharing_hist_Y_vec[i + 1] + sharing_v_Y_vec[i + 1]) & mask_hist;
    }

    std::vector<ShareA> sharing_w_X_vec(sharing_hist_X_vec);
    std::vector<ShareA> sharing_w_Y_vec(sharing_hist_Y_vec);
    std::transform(std::execution::par_unseq,
                   sharing_w_X_vec.begin(), sharing_w_X_vec.end(), sharing_v_X_vec.begin(),
                   sharing_w_X_vec.begin(), [mask_hist](auto wx, auto vx)
                   { return (wx + vx * 2) & mask_hist; });
    std::transform(std::execution::par_unseq,
                   sharing_w_Y_vec.begin(), sharing_w_Y_vec.end(), sharing_v_Y_vec.begin(),
                   sharing_w_Y_vec.begin(), [mask_hist](auto wy, auto vy)
                   { return (wy + vy * 2) & mask_hist; });

    std::vector<ShareA> sharing_hist_YX_vec = my::append_vectors(sharing_hist_Y_vec, sharing_hist_X_vec);
    std::vector<ShareA> sharing_w_XY_vec = my::append_vectors(sharing_w_X_vec, sharing_w_Y_vec);
    std::vector<ShareA> sharing_u_XY_vec(2 * len_hist);

    // 多线程乘法
    std::thread secure_hadamard[num_threads];
    size_t chunk_size = 2 * len_hist / num_threads;
    for (int tid = 0; tid < num_threads; ++tid)
    {
        int offset = tid * chunk_size;
        int num_ops;
        if (tid == (num_threads - 1))
        {
            num_ops = 2 * len_hist - offset;
        }
        else
        {
            num_ops = chunk_size;
        }
        // 结果放大了 scale*2*2,
        // if(bw_hist<required_bit_width){
        //     secure_hadamard[tid] = std::thread(
        //         my::hadamard_product_thread, party, iopackArr, otpackArr, tid,
        //         num_ops, sharing_hist_YX_vec.data() + offset, sharing_w_XY_vec.data() + offset, sharing_u_XY_vec.data() + offset,
        //         bw_hist, bw_hist, bw_numeric, false, false, MultMode::None, nullptr, nullptr);
        // }
        // else{
        secure_hadamard[tid] = std::thread(
            my::hadamard_product_thread, party, iopackArr, otpackArr, tid,
            num_ops, sharing_hist_YX_vec.data() + offset, sharing_w_XY_vec.data() + offset, sharing_u_XY_vec.data() + offset,
            bw_hist, bw_hist, bw_hist, false, false, MultMode::None, nullptr, nullptr);
        // }
    }
    for (int tid = 0; tid < num_threads; ++tid)
    {
        secure_hadamard[tid].join();
    }

#ifdef _PERFORMANCE_
    uint64_t hadamard_product_rounds = iopackArr[0]->get_rounds() - last_rounds;
    last_rounds = iopackArr[0]->get_rounds();
    fmt::print(fg(fmt::color::white_smoke), "{}->hadamard_product_rounds= {}, num rounds= {}\n", party_name, hadamard_product_rounds, last_rounds);
    logger->info("{}->hadamard_product_rounds= {}, num rounds= {}", party_name, hadamard_product_rounds, last_rounds);
#endif

    ShareA sharing_U_X = std::reduce(std::execution::par_unseq, sharing_u_XY_vec.begin(), sharing_u_XY_vec.begin() + len_hist) & mask_hist;
    ShareA sharing_U_Y = std::reduce(std::execution::par_unseq, sharing_u_XY_vec.begin() + len_hist, sharing_u_XY_vec.end()) & mask_hist;
    std::vector<ShareA> sharing_U_XY_vec = {sharing_U_X, sharing_U_Y};

    /**** Process & Write Benchmarking Data *****/
#ifdef _PERFORMANCE_
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
    auto hist_X_vec = my::reconstruct_uint64_t(party, sci::ALICE, iopackArr, sharing_hist_X_vec, bw_hist, __LINE__);
    auto hist_Y_vec = my::reconstruct_uint64_t(party, sci::ALICE, iopackArr, sharing_hist_Y_vec, bw_hist, __LINE__);
    auto hist_YX_vec = my::reconstruct_uint64_t(party, sci::ALICE, iopackArr, sharing_hist_YX_vec, bw_hist, __LINE__);
    auto w_X_vec = my::reconstruct_uint64_t(party, sci::ALICE, iopackArr, sharing_w_X_vec, bw_hist, __LINE__);
    auto w_Y_vec = my::reconstruct_uint64_t(party, sci::ALICE, iopackArr, sharing_w_Y_vec, bw_hist, __LINE__);
    auto w_XY_vec = my::reconstruct_uint64_t(party, sci::ALICE, iopackArr, sharing_w_XY_vec, bw_hist, __LINE__);
    auto u_XY_vec = my::reconstruct_uint64_t(party, sci::ALICE, iopackArr, sharing_u_XY_vec, bw_hist, __LINE__);
#endif

    auto U_XY_vec = my::reconstruct_uint64_t(party, sci::ALICE, iopackArr, sharing_U_XY_vec, bw_hist);

    if (party == ALICE)
    {
        fmt::print(fg(fmt::color::blue), "-----------------------Verification in SecChiQuared-----------------------\n");
        fmt::print(fg(fmt::color::blue), "{}->[U_X, U_Y]= {}\n", party_name, U_XY_vec);
        fmt::print(fg(fmt::color::blue), "{}->bw_hist= {}\n", party_name, bw_hist);
#ifdef _CHECK_
        fmt::print("{}->hist_X_vec= {}\n", party_name, hist_X_vec);
        fmt::print("{}->hist_Y_vec= {}\n", party_name, hist_Y_vec);
        fmt::print("{}->hist_YX_vec= {}\n", party_name, hist_YX_vec);
        fmt::print("{}->w_X_vec= {}\n", party_name, w_X_vec);
        fmt::print("{}->w_Y_vec= {}\n", party_name, w_Y_vec);
        fmt::print("{}->w_XY_vec= {}\n", party_name, w_XY_vec);
        fmt::print("{}->u_XY_vec= {}\n", party_name, u_XY_vec);
#endif
    }
#endif

    /******************* Cleanup ****************/
    spdlog::drop_all();
    return sharing_U_XY_vec;
}

int main(int argc, char **argv)
{
    int left = -1;
    int step = -1;
    int number_of_data = -1;
    uint64_t generate_min = -1;
    uint64_t generate_max = -1;
    int len_hist = -1;
    int bw_hist = -1;
    int bw_numeric = -1;

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
    amap.parse(argc, argv);

    my::set_value_int(number_of_data, 500);
    my::set_value_int(len_hist, 100);
    my::set_value_int(left, 100);
    my::set_value_int(step, 1);
    my::set_value_uint64_t(generate_min, left + len_hist * step * 0.02);
    my::set_value_uint64_t(generate_max, left + (len_hist - 1) * step - len_hist * step * 0.02);
    
    my::set_value_int(bw_hist, ceil(2*log2(number_of_data)+log2(3*len_hist)));
    
    // my::set_value_int(bw_hist, 2 * (ceil(log2(number_of_data)) + 2));
    // my::set_value_int(bw_hist, (ceil(log2(number_of_data)) + 2));
    my::set_value_int(bw_numeric, ceil(log2(left + (len_hist - 1) * step)));

    log_filename = fmt::format("logs/my_ring_U_{}.txt",party);

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
               "-------------- {}-> U ----------------------\nnumber_of_data:{}, bw_numeric:{}, len_hist:{}, bw_hist:{}, left:{}, generate_min:{}, generate_max:{}\n",
               party_name, number_of_data, bw_numeric, len_hist, bw_hist, left, generate_min, generate_max);
    logger->info("-------------- {}-> U ----------------------\nnumber_of_data:{}, bw_numeric:{}, len_hist:{}, bw_hist:{}, left:{}, generate_min:{}, generate_max:{}\n",
                 party_name, number_of_data, bw_numeric, len_hist, bw_hist, left, generate_min, generate_max);

    my::setup_io_and_ots(party, iopackArr, otpackArr, num_threads, port, address);
    fmt::print("All Base OTs Done\n");

    /************ Generate Test Data ************/

    std::vector<uint64_t> hist_X_vec_actual_ALICE;
    std::vector<uint64_t> hist_Y_vec_actual_ALICE;
    std::vector<uint64_t> sharing_hist_X_vec;
    std::vector<uint64_t> sharing_hist_Y_vec;

    if (party == ALICE)
    {
        // std::vector<uint64_t> numeric_data_X_vec{5, 6, 11, 10, 12, 46, 42};
        // std::vector<uint64_t> numeric_data_Y_vec{50, 61, 13, 9, 74, 45, 75};

        auto numeric_data_X_vec = my::generate_numeric_data_vector(number_of_data, bw_numeric, generate_min, generate_max);
        auto numeric_data_Y_vec = my::generate_numeric_data_vector(number_of_data, bw_numeric, generate_min, generate_max);

        hist_X_vec_actual_ALICE = my::numeric_data_to_hist(len_hist, left, step, numeric_data_X_vec, number_of_data);
        hist_Y_vec_actual_ALICE = my::numeric_data_to_hist(len_hist, left, step, numeric_data_Y_vec, number_of_data);

        sharing_hist_X_vec = my::generate_and_send_sharing_vector_plaintext(iopackArr[0], hist_X_vec_actual_ALICE, bw_hist);
        sharing_hist_Y_vec = my::generate_and_send_sharing_vector_plaintext(iopackArr[0], hist_Y_vec_actual_ALICE, bw_hist);
    }
    else if (party == BOB)
    {
        sharing_hist_X_vec = my::recv_sharing_vector_plaintext_vector(iopackArr[0], len_hist);
        sharing_hist_Y_vec = my::recv_sharing_vector_plaintext_vector(iopackArr[0], len_hist);
    }

    auto sharing_U_XY_vec = Secure_U_Test(sharing_hist_X_vec, sharing_hist_Y_vec, len_hist, bw_hist);

    /************** Verification ****************/

#ifdef _CHECK_
    if (party == BOB)
    {
        // iopackArr[0]->io->send_data(&sharing_Q_1_4, sizeof(ShareA));
    }
    else if (party == ALICE)
    {
        fmt::print("-----------------------Verification-----------------------\n");

        // fmt::print("{}->hist_vec_actual_ALICE: {}\n", party_name, hist_vec_actual_ALICE);
        // fmt::print("{}->R_mat_M: {}\n", party_name, R_mat_M);
        // fmt::print("{}->P_mat_L: {}\n", party_name, P_mat_L);
    }
#endif

    for (int i = 0; i < num_threads; i++)
    {
        delete iopackArr[i];
        delete otpackArr[i];
    }
    spdlog::drop_all();
}
