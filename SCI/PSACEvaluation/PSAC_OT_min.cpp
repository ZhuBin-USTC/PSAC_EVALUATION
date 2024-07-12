// We are just using multiple rounds of comparison operations to model the overhead of the sorted network algorithm, not really implementing it.

#include "Math/math-functions.h"
#include "MyFunctions/my-functions.h"
#include <fstream>
#include <iostream>
#include <thread>
#include <fmt/format.h>
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

uint64_t hist_min(
    std::vector<ShareA> &sharing_hist_vec,
    int bw_hist,
    int len_hist,
    uint64_t data_space_left, // 明文空间左边界
    uint64_t data_space_step, // 明文空间增量
    int bw_numeric            // 返回值的bit长度，也是明文空间的bit长度
)
{
    assert(bw_numeric >= ceil(log2(data_space_left + (len_hist - 1) * data_space_step)));

#ifdef _PERFORMANCE_

    std::shared_ptr<spdlog::logger> logger;
    try
    {
        logger = spdlog::basic_logger_st<spdlog::async_factory>("min", log_filename);
    }
    catch (const spdlog::spdlog_ex &ex)
    {
        std::cout << "Log init failed: " << ex.what() << std::endl;
    }
    string party_name = (party == ALICE ? "ALICE" : "BOB");
    auto rounds_begin = iopackArr[0]->get_rounds();
    fmt::print(fg(fmt::color::green), "{}->rounds_begin= {}, bw_hist= {}, bw_compare= {}\n", party_name, rounds_begin,bw_hist,bw_numeric);
    logger->info("{}->rounds_begin= {}, bw_hist= {}, bw_compare= {}", party_name, rounds_begin,bw_hist,bw_numeric);

    uint64_t total_comm = 0;
    uint64_t thread_comm[num_threads];
    for (int i = 0; i < num_threads; i++)
    {
        thread_comm[i] = iopackArr[i]->get_comm();
    }
    auto start = sci::clock_start();
#endif
    /******************************************************/

    uint64_t mask_hist = my::mask_gen_uint64_t(bw_hist);
    uint64_t mask_numeric = my::mask_gen_uint64_t(bw_numeric);

    // TODO 修改1
    auto sharing_v_vec = std::vector<ShareA>(len_hist);
    sharing_v_vec[0] = sharing_hist_vec[0];
    for (int i = 1; i < len_hist; i++)
    {
        sharing_v_vec[i] = (sharing_v_vec[i - 1] + sharing_hist_vec[i]) & mask_hist;
    }

    auto sharing_eq_vec_B = std::vector<ShareB>(len_hist); // 比较大小返回值向量
    auto sharing_eq_vec_A = std::vector<ShareA>(len_hist);

    /************** Fork Threads ****************/
    std::thread compare_threads_arithmetic[num_threads];
    int chunk_size = len_hist / num_threads;
    for (int i = 0; i < num_threads; ++i)
    {
        int offset = i * chunk_size;
        int num_ops;
        if (i == (num_threads - 1))
        {
            num_ops = len_hist - offset;
        }
        else
        {
            num_ops = chunk_size;
        }

        const uint64_t compare_value = 0;
        compare_threads_arithmetic[i] = std::thread(
            my::compare_arithmetic_thread, EQUAL, party, iopackArr, otpackArr, i,
            sharing_v_vec.data() + offset, compare_value, num_ops, bw_hist, bw_numeric, nullptr,
            sharing_eq_vec_B.data() + offset, nullptr, sharing_eq_vec_A.data() + offset);
    }
    for (int i = 0; i < num_threads; ++i)
    {
        compare_threads_arithmetic[i].join();
    }

    uint64_t compare_rounds = iopackArr[0]->get_rounds() - rounds_begin;

    uint64_t sharing_count_equal = std::reduce(std::execution::par_unseq, sharing_eq_vec_A.begin(), sharing_eq_vec_A.end(), 0, std::plus<ShareA>()) & mask_numeric;

    uint64_t sharing_min_numeric = sharing_count_equal * data_space_step;
    sharing_min_numeric = (party == ALICE ? sharing_min_numeric + data_space_left : sharing_min_numeric);
    sharing_min_numeric = sharing_min_numeric & mask_numeric;

#ifdef _PERFORMANCE_

    long long t = time_from(start);
    for (int i = 0; i < num_threads; i++)
    {
        thread_comm[i] = iopackArr[i]->get_comm() - thread_comm[i];
        total_comm += thread_comm[i];
    }
    /**** Process & Write Benchmarking Data *****/
    fmt::print(fg(fmt::color::white_smoke), "{}->compare rounds= {}, num rounds= {}\n", party_name, compare_rounds, iopackArr[0]->get_rounds());
    logger->info("{}->compare rounds= {}, num rounds= {}", party_name, compare_rounds, iopackArr[0]->get_rounds());

    fmt::print(fg(fmt::color::white_smoke), "{}->Number of dim/s:\t{}\n", party_name, (double(len_hist) / t) * 1e6);
    fmt::print(fg(fmt::color::white_smoke), "{}->Computing Time\t{} ms\n", party_name, t / (1000.0));
    fmt::print(fg(fmt::color::white_smoke), "{}->Computing Bytes Sent\t{} bytes, {} MB\n", party_name, total_comm, total_comm / (1024.0 * 1024.0));
    logger->info("{}->Number of dim/s:\t{}", party_name, (double(len_hist) / t) * 1e6);
    logger->info("{}->Computing Time\t{} ms", party_name, t / (1000.0));
    logger->info("{}->Computing Bytes Sent\t{} bytes, {} MB", party_name, total_comm, total_comm / (1024.0 * 1024.0));
#endif
    spdlog::drop_all();
    return sharing_min_numeric;
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

    my::set_value_int(len_hist, 10000);
    my::set_value_int(left, 10);
    my::set_value_int(step, 1);

    my::set_value_int(number_of_data, 1000);
    my::set_value_uint64_t(generate_min, left + len_hist * step * 0.02);
    my::set_value_uint64_t(generate_max, left + (len_hist - 1) * step - len_hist * step * 0.02);

    // my::set_value_int(bw_hist, 16);
    my::set_value_int(bw_hist, ceil(log2(number_of_data)));
    my::set_value_int(bw_numeric, ceil(log2(left + (len_hist - 1) * step)));

    log_filename = fmt::format("logs/my_ring_min_{}.txt",party);

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
               "-------------- {}-> min ----------------------\nnumber_of_data:{}, bw_numeric:{}, len_hist:{}, bw_hist:{}, left:{}, generate_min:{}, generate_max:{}\n",
               party_name, number_of_data, bw_numeric, len_hist, bw_hist, left, generate_min, generate_max);
    logger->info("-------------- {}-> min ----------------------\nnumber_of_data:{}, bw_numeric:{}, len_hist:{}, bw_hist:{}, left:{}, generate_min:{}, generate_max:{}\n",
               party_name, number_of_data, bw_numeric, len_hist, bw_hist, left, generate_min, generate_max);

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
    std::vector<uint64_t> numeric_data_vec;
    std::vector<uint64_t> sharing_hist_vec;

    if (party == ALICE)
    {
        numeric_data_vec = my::generate_numeric_data_vector(number_of_data, bw_numeric, generate_min, generate_max);
        auto hist_vec = my::numeric_data_to_hist(len_hist, left, step, numeric_data_vec, number_of_data);
        sharing_hist_vec = my::generate_and_send_sharing_vector_plaintext(iopackArr[0], hist_vec, bw_hist);
    }
    else if (party == BOB)
    {
        sharing_hist_vec = my::recv_sharing_vector_plaintext_vector(iopackArr[0], len_hist);
    }

    uint64_t sharing_min_numeric = hist_min(sharing_hist_vec, bw_hist, len_hist, left, step, bw_numeric);

/************** Verification ****************/
#ifdef _VERIFY_
    std::vector<ShareA> sharing_min_numeric_vec(1, sharing_min_numeric);
    auto min_numeric_vec = my::reconstruct_uint64_t(party, sci::ALICE, iopackArr, sharing_min_numeric_vec, bw_numeric);
    if (party == ALICE)
    {
        uint64_t min_actual = *(std::min_element(numeric_data_vec.data(), numeric_data_vec.data() + number_of_data));
        fmt::print(fg(fmt::color::blue), "Line:{}, min={} min_actual: {}\n", __LINE__, min_numeric_vec, min_actual);
        assert((min_actual - left) / step == (min_numeric_vec[0] - left) / step);
    }
    else if (party == BOB)
    {
        iopackArr[0]->io->send_data(&sharing_min_numeric, sizeof(uint64_t));
    }
#endif

    for (int i = 0; i < num_threads; i++)
    {
        delete iopackArr[i];
        delete otpackArr[i];
    }
    spdlog::drop_all();
}
