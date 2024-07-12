#include "Math/math-functions.h"
#include <iostream>
#include <thread>
#include "MyFunctions/my-functions.h"
#include <vector>
#include <algorithm>

#include <fmt/core.h>
#include <fmt/color.h>

#include "library_fixed.h"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/async.h"

using namespace sci;

#define MAX_THREADS 4

int party, port = 32000;
int num_threads = 1;
string address = "127.0.0.1";
spdlog::filename_t log_filename;

int32_t bitlength = 32;

// test 01
uint64_t ProtocolTimeInMilliSec = 0;
uint64_t ProtocolCommSent = 0;

void rmind_min(
    int num_data,
    std::unique_ptr<uint64_t[]> &sharing_data_ptr,
    // std::unique_ptr<uint64_t[]> &res_data_ptr,
    int bitlength_arithmetic)
{

// test 02
#ifdef LOG_LAYERWISE
    INIT_TIMER;
    INIT_ALL_IO_DATA_SENT;
#endif

    std::shared_ptr<spdlog::logger> logger;
    try
    {
        logger = spdlog::basic_logger_st<spdlog::async_factory>("min", log_filename);
    }
    catch (const spdlog::spdlog_ex &ex)
    {
        std::cout << "Log init failed: " << ex.what() << std::endl;
    }

    fmt::print(fg(fmt::color::red), "Begin Number of rounds = {}\n", iopack->get_rounds() - num_rounds);
    logger->info("Begin Number of rounds = {}", iopack->get_rounds() - num_rounds);
    
    auto res_data_ptr = std::make_unique<uint64_t[]>(num_data);
    my::VectorSortEmulate(num_data, sharing_data_ptr.get(), res_data_ptr.get(), bitlength_arithmetic);

    fmt::print(fg(fmt::color::red), "VectorSortEmulate Number of rounds = {}\n", iopack->get_rounds() - num_rounds);
    logger->info("VectorSortEmulate Number of rounds = {}", iopack->get_rounds() - num_rounds);

    // int K = ceil(num_data / 256.0);
    // uint64_t **spec = new uint64_t *[K];

    // for (int i = 0; i < K - 1; i++)
    // {
    //     spec[i] = new uint64_t[256];
    //     memcpy(spec[i], sharing_data_ptr.get() + i * 256, 256 * sizeof(uint64_t));
    // }

    // spec[K - 1] = new uint64_t[256];
    // memcpy(spec[K - 1], sharing_data_ptr.get() + (K - 1) * 256, (num_data - (K - 1) * 256) * sizeof(uint64_t));

    // for (int i = num_data - (K - 1) * 256; i < 256; i++)
    // {
    //     spec[K - 1][i] = 0;
    // }

    // int KK = ceil(log2(K));
    // int KK_power2 = pow(2, KK);
    // uint64_t *x = new uint64_t[KK_power2];
    // for (int i = 0; i < KK_power2; i++)
    // {
    //     x[i] = 10;
    // }
    // uint64_t *y = new uint64_t[KK_power2];
    // my::MultipleLUT(spec, x, y, K, 8, bitlength_arithmetic);

    // if (K > 1)
    // {
    //     uint64_t **spec2 = new uint64_t *[1];
    //     spec2[0] = y;
    //     uint64_t x2 = 0;
    //     uint64_t y2 = 0;
    //     my::MultipleLUT(spec2, &x2, &y2, 1, KK, bitlength_arithmetic);
    //     delete[] spec2;
    // }

// #ifdef LOG_LAYERWISE
//     fmt::print(fg(fmt::color::red), "MultipleLUT Number of rounds = {}\n", iopack->get_rounds() - num_rounds);
//     logger->info("MultipleLUT Number of rounds = {}", iopack->get_rounds() - num_rounds);
// #endif

//     delete[] x;
//     delete[] y;
//     for (int i = 0; i < K; i++)
//     {
//         delete[] spec[i];
//     }
//     delete[] spec;

// test 03
#ifdef LOG_LAYERWISE
    auto temp = TIMER_TILL_NOW;
    ProtocolTimeInMilliSec = temp;
    uint64_t curComm;
    FIND_ALL_IO_TILL_NOW(curComm);
    ProtocolCommSent = curComm;

    fmt::print(fg(fmt::color::red), "SortEmulateTimeInMilliSec= {}, SortEmulateCommSent= {} MiB\n", SortEmulateTimeInMilliSec, ((SortEmulateCommSent) / (1.0 * (1ULL << 20))));
    // fmt::print(fg(fmt::color::red), "MultipleLUTTimeInMilliSec= {}, MultipleLUTCommSent= {} MiB\n", MultipleLUTTimeInMilliSec, ((MultipleLUTCommSent) / (1.0 * (1ULL << 20))));
    fmt::print(fg(fmt::color::red), "ProtocolTimeInMilliSec= {}, ProtocolCommSent= {} MiB\n", ProtocolTimeInMilliSec, ((ProtocolCommSent) / (1.0 * (1ULL << 20))));

    logger->info("SortEmulateTimeInMilliSec= {}, SortEmulateCommSent= {} MiB", SortEmulateTimeInMilliSec, ((SortEmulateCommSent) / (1.0 * (1ULL << 20))));
    // logger->info("MultipleLUTTimeInMilliSec= {}, MultipleLUTCommSent= {} MiB", MultipleLUTTimeInMilliSec, ((MultipleLUTCommSent) / (1.0 * (1ULL << 20))));
    logger->info("ProtocolTimeInMilliSec= {}, ProtocolCommSent= {} MiB", ProtocolTimeInMilliSec, ((ProtocolCommSent) / (1.0 * (1ULL << 20))));

    if (party == SERVER)
    {
        uint64_t ProtocolCommSentClient = 0;
        io->recv_data(&ProtocolCommSentClient, sizeof(uint64_t));
        fmt::print(
            fg(fmt::color::red), 
            "ProtocolComm(sent+received) = {} MiB\n", 
            ((ProtocolCommSent + ProtocolCommSentClient) / (1.0 * (1ULL << 20))));
        logger->info(
            "ProtocolComm(sent+received) = {} MiB", 
            ((ProtocolCommSent + ProtocolCommSentClient) / (1.0 * (1ULL << 20))));
    }
    else if (party == CLIENT)
    {
        io->send_data(&ProtocolCommSent, sizeof(uint64_t));
    }
#endif
}

int main(int argc, char **argv)
{
    int len_hist = -1;
    int left = -1;
    int step = -1;
    int number_of_data = -1;
    uint64_t generate_min = -1;
    uint64_t generate_max = -1;
    int bw_numeric = -1;

    /************* Argument Parsing  ************/
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.arg("po", port, "Port Number");
    amap.arg("nt", num_threads, "Number of threads");
    amap.arg("lh", len_hist, "Length of histgram vector");
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

    my::set_value_int(bw_numeric, ceil(log2(left + (len_hist - 1) * step)) + 1);
    bitlength = bw_numeric;

    dim_frequency_global = len_hist;
    num_data_global = number_of_data;

    fmt::print(fg(fmt::color::tomato), "generate_min={},generate_max={}\n", generate_min, generate_max);
    fmt::print(fg(fmt::color::tomato), "bw_arithmetic={}\n", bitlength);

    log_filename = fmt::format("logs/my_ring_min_rmind_OT_{}.txt", party);

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
    fmt::print(
        fg(fmt::color::blue),
        "-------------- {}-> min rmind---------------\nnumber_of_data:{}, bw_numeric:{}, len_hist:{}, left:{}, generate_min:{}, generate_max:{}\n",
        party_name, number_of_data, bw_numeric, len_hist, left, generate_min, generate_max);
    logger->info(
        "-------------- {}-> min rmind---------------\nnumber_of_data:{}, bw_numeric:{}, len_hist:{},  left:{}, generate_min:{}, generate_max:{}\n",
        party_name, number_of_data, bw_numeric, len_hist, left, generate_min, generate_max);

    StartComputation();
    std::unique_ptr<uint64_t[]> sharing_data_ptr;

    auto iopackMain = new IOPack(party, port - 1, address);

    if (party == ALICE)
    {
        auto data_arithmetic_ptr = my::generated_worker_data_arithmetic(number_of_data, bitlength, generate_min, generate_max);
        sharing_data_ptr = my::generate_and_send_sharing_vector_plaintext(iopackMain, data_arithmetic_ptr, number_of_data, bitlength);
    }
    else if (party == BOB)
    {
        sharing_data_ptr = my::recv_sharing_vector_plaintext(iopackMain, number_of_data);
    }
    delete iopackMain;

    rmind_min(number_of_data, sharing_data_ptr, bitlength);
    EndComputation();


    spdlog::drop_all();
}
