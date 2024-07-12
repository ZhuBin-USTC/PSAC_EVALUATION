// We use the bitonic_sort function provided by SCI to evaluate the performance of the sorting algorithm using GC.

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
#include "GC/emp-sh2pc.h"
// #include <memory>
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

int size_1024 = 64;

void sort_gc(
    int num_data,
    std::unique_ptr<uint64_t[]> &sharing_data_ptr,
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

    setup_semi_honest(iopackArr[0]->io_GC, party, 1024 * size_1024);

    auto tmp = iopackArr[0]->get_comm();
    auto tmp2 = iopackArr[0]->get_rounds();

    auto zero = Integer(bitlength_arithmetic, 0, PUBLIC);

    auto shareA = std::make_unique<Integer[]>(num_data);
    auto shareB = std::make_unique<Integer[]>(num_data);
    auto Input = std::make_unique<Integer[]>(num_data);

    sci::PRG128 prg;
    auto share_res_A = std::make_unique<uint64_t[]>(num_data);
    auto share_res_B = std::make_unique<uint64_t[]>(num_data);
    auto share_res_A_Integer = std::make_unique<Integer[]>(num_data);
    auto share_res_B_Integer = std::make_unique<Integer[]>(num_data);

    uint64_t mask_data = (bitlength_arithmetic == 64 ? -1 : ((1ULL << bitlength_arithmetic) - 1));
    prg.random_data(share_res_A.get(), num_data * sizeof(uint64_t));
    for (int i = 0; i < num_data; i++)
    {
        share_res_A[i] = share_res_A[i] & mask_data;
    }
    for (int i = 0; i < num_data; i++)
    {
        share_res_A_Integer[i] = Integer(bitlength_arithmetic, share_res_A[i], ALICE);
    }

    for (int i = 0; i < num_data; i++)
    {
        shareA[i] = Integer(bitlength_arithmetic, sharing_data_ptr[i], ALICE);
    }
    for (int i = 0; i < num_data; i++)
    {
        shareB[i] = Integer(bitlength_arithmetic, sharing_data_ptr[i], BOB);
    }

    for (int i = 0; i < num_data; i++)
    {
        Input[i] = shareA[i] + shareB[i];
    }

    sort<Integer>(Input.get(), num_data);

    fmt::print(
        fg(fmt::color::yellow_green),
        "sort dim_frequency={}, GC: {} MiB, {} rounds\n",
        num_data,
        (iopackArr[0]->get_comm() - tmp) / (1.0 * (1ULL << 20)), iopackArr[0]->get_rounds() - tmp2);

    logger->info(
        "sort dim_frequency={}, GC: {} MiB, {} rounds",
        num_data,
        (iopackArr[0]->get_comm() - tmp) / (1.0 * (1ULL << 20)), iopackArr[0]->get_rounds() - tmp2);

    for (int i = 0; i < num_data; i++)
    {
        share_res_B_Integer[i] = Input[i] - share_res_A_Integer[i];
    }

    for (int i = 0; i < num_data; i++)
    {
        share_res_B[i] = share_res_B_Integer[i].reveal<uint64_t>(BOB);
    }

// test 03
#ifdef LOG_LAYERWISE
    auto temp = TIMER_TILL_NOW;
    ProtocolTimeInMilliSec = temp;
    uint64_t curComm;
    FIND_ALL_IO_TILL_NOW(curComm);
    ProtocolCommSent = curComm;

    fmt::print(
        fg(fmt::color::green),
        "Sharing results dim_frequency={}, GC: {} MiB, {} rounds\n",
        num_data, (iopackArr[0]->get_comm() - tmp) / (1.0 * (1ULL << 20)), iopackArr[0]->get_rounds() - tmp2);
    fmt::print(fg(fmt::color::red), "AND Gates = {}\n", circ_exec->num_and());
    fmt::print(fg(fmt::color::green), "Time in ms for current Protocol = {}, ProtocolCommSent={} MiB\n", ProtocolTimeInMilliSec, ((ProtocolCommSent) / (1.0 * (1ULL << 20))));

    logger->info(
        "Sharing results dim_frequency={}, GC: {} MiB, {} rounds",
        num_data, (iopackArr[0]->get_comm() - tmp) / (1.0 * (1ULL << 20)), iopackArr[0]->get_rounds() - tmp2);
    logger->info("AND Gates = {}", circ_exec->num_and());
    logger->info("Time in ms for current Protocol = {}, ProtocolCommSent={} MiB", ProtocolTimeInMilliSec, ((ProtocolCommSent) / (1.0 * (1ULL << 20))));
#endif

    Bit x(false, PUBLIC);
    fmt::print("finish.\n", x.reveal());

    // fmt::print(
    //     fg(fmt::color::blue),
    //     "ProtocolTimeInMilliSec={}, ProtocolCommSent={} MiB\n", ProtocolTimeInMilliSec, ((ProtocolCommSent) / (1.0 * (1ULL << 20))));
    // logger->info(
    //     "ProtocolTimeInMilliSec={}, ProtocolCommSent={} MiB",
    //     ProtocolTimeInMilliSec, ((ProtocolCommSent) / (1.0 * (1ULL << 20))));

    // if (party == SERVER)
    // {
    //     uint64_t ProtocolCommSentClient = 0;
    //     io->recv_data(&ProtocolCommSentClient, sizeof(uint64_t));
    //     fmt::print(
    //         fg(fmt::color::blue_violet),
    //         "ProtocolComm(sent+received) = {} MiB\n",
    //         ((ProtocolCommSent + ProtocolCommSentClient) / (1.0 * (1ULL << 20))));
    //     logger->info(
    //         "ProtocolComm(sent+received) = {} MiB",
    //         ((ProtocolCommSent + ProtocolCommSentClient) / (1.0 * (1ULL << 20))));
    // }
    // else if (party == CLIENT)
    // {
    //     io->send_data(&ProtocolCommSent, sizeof(uint64_t));
    // }
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
    amap.arg("sgc", size_1024, "gc batch size");
    amap.parse(argc, argv);

    my::set_value_int(len_hist, 10000);
    my::set_value_int(left, 10);
    my::set_value_int(step, 1);

    my::set_value_int(number_of_data, 1000);
    my::set_value_uint64_t(generate_min, left + len_hist * step * 0.02);
    my::set_value_uint64_t(generate_max, left + (len_hist - 1) * step - len_hist * step * 0.02);

    // l''=2, p'=1
    my::set_value_int(bw_numeric, ceil(log2(left + (len_hist - 1) * step)) + 2);
    bitlength = bw_numeric;

    fmt::print(fg(fmt::color::tomato), "generate_min={},generate_max={}\n", generate_min, generate_max);
    fmt::print(fg(fmt::color::tomato), "bw_arithmetic={}\n", bitlength);

    log_filename = fmt::format("logs/my_gc_min_rmind_{}.txt", party);

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
        "-------------- {}-> min rmind gc----------\nnumber_of_data:{}, bw_numeric:{}, len_hist:{}, left:{}, generate_min:{}, generate_max:{}, sgc={}\n",
        party_name, number_of_data, bw_numeric, len_hist, left, generate_min, generate_max,size_1024);
    logger->info(
        "-------------- {}-> min rmind gc----------\nnumber_of_data:{}, bw_numeric:{}, len_hist:{},  left:{}, generate_min:{}, generate_max:{}, sgc={}\n",
        party_name, number_of_data, bw_numeric, len_hist, left, generate_min, generate_max,size_1024);

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

    sort_gc(number_of_data, sharing_data_ptr, bitlength);
    // EndComputation();


    spdlog::drop_all();
}
