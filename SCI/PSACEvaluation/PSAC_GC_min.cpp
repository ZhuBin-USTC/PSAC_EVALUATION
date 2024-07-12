#include "Math/math-functions.h"
#include <iostream>
#include <thread>
#include "MyFunctions/my-functions.h"
#include <vector>
#include <algorithm>

#include <fmt/core.h>
#include <fmt/color.h>

#include "library_fixed.h"
#include "GC/emp-sh2pc.h"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/async.h"

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

int size_1024 = 16;

void hist_min2(
    std::vector<ShareA> &sharing_hist_vec,
    int bw_hist,
    int len_hist,
    uint64_t data_space_left, // 明文空间左边界
    uint64_t data_space_step, // 明文空间增量
    int bw_numeric            // 返回值的bit长度，也是明文空间的bit长度
)
{
#ifdef LOG_LAYERWISE
    INIT_TIMER;
    INIT_ALL_IO_DATA_SENT;
#endif

    uint64_t mask_hist = my::mask_gen_uint64_t(bw_hist);
    auto sharing_v_vec = std::vector<ShareA>(len_hist);
    sharing_v_vec[0] = sharing_hist_vec[0];
    for (int i = 1; i < len_hist; i++)
    {
        sharing_v_vec[i] = (sharing_v_vec[i - 1] + sharing_hist_vec[i]) & mask_hist;
    }
    std::cout << std::endl;
    setup_semi_honest(iopackArr[0]->io_GC, party, 1024 * size_1024);

    auto tmp = iopackArr[0]->get_comm();
    auto tmp2 = iopackArr[0]->get_rounds();

    auto shareA = std::make_unique<Integer[]>(len_hist);
    auto shareB = std::make_unique<Integer[]>(len_hist);
    auto compare_res_bool = std::make_unique<Bit[]>(len_hist);
    // auto compare_res_Integer = std::make_unique<Integer[]>(len_hist);
    auto compare_res_share_A_Integer = std::make_unique<Integer[]>(len_hist);
    auto compare_res_share_B_Integer = std::make_unique<Integer[]>(len_hist);

    sci::PRG128 prg;
    auto share_res_A = std::make_unique<uint64_t[]>(len_hist);
    auto share_res_B = std::make_unique<uint64_t[]>(len_hist);
    prg.random_data(share_res_A.get(), len_hist * sizeof(uint64_t));
    for (int i = 0; i < len_hist; i++)
    {
        share_res_A[i] = share_res_A[i] & mask_hist;
    }

    auto zero = Integer(bw_hist, 0, PUBLIC);
    for (int i = 0; i < len_hist; i++)
    {
        shareA[i] = Integer(bw_hist, sharing_hist_vec[i], ALICE);
    }
    for (int i = 0; i < len_hist; i++)
    {
        shareB[i] = Integer(bw_hist, sharing_hist_vec[i], BOB);
    }
    for (int i = 0; i < len_hist; i++)
    {
        compare_res_bool[i] = ((shareA[i] + shareB[i]) == zero);
    }

    fmt::print(fg(fmt::color::red), "Equal->Bit len_hist={}, GC: {} MiB, {} rounds\n",
               len_hist,
               (iopackArr[0]->get_comm() - tmp) / (1.0 * (1ULL << 20)),
               iopackArr[0]->get_rounds() - tmp2);

    for (int i = 0; i < len_hist; i++)
    {
        compare_res_share_B_Integer[i] = If(
            compare_res_bool[i],
            Integer(bw_hist, (1 - share_res_A[i]) & mask_hist, ALICE),
            Integer(bw_hist, (-share_res_A[i]) & mask_hist, ALICE));
    }

    fmt::print(fg(fmt::color::red), "Equal->Integer len_hist={}, GC: {} MiB, {} rounds\n",
               len_hist,
               (iopackArr[0]->get_comm() - tmp) / (1.0 * (1ULL << 20)),
               iopackArr[0]->get_rounds() - tmp2);

    for (int i = 0; i < len_hist; i++)
    {
        share_res_B[i] = compare_res_share_B_Integer[i].reveal<uint64_t>(BOB);
    }

    fmt::print(fg(fmt::color::red), "Equal->Share len_hist={}, GC: {} MiB, {} rounds\n",
               len_hist,
               (iopackArr[0]->get_comm() - tmp) / (1.0 * (1ULL << 20)),
               iopackArr[0]->get_rounds() - tmp2);
    fmt::print(fg(fmt::color::red), "AND Gates = {}\n", circ_exec->num_and());

// test 03
#ifdef LOG_LAYERWISE
    auto temp = TIMER_TILL_NOW;
    ProtocolTimeInMilliSec = temp;
    fmt::print(fg(fmt::color::red), "Time in sec for current Protocol = {}, ProtocolCommSent={}\n",
               (ProtocolTimeInMilliSec / 1000.0),
               ((ProtocolCommSent) / (1.0 * (1ULL << 20))));
    uint64_t curComm;
    FIND_ALL_IO_TILL_NOW(curComm);
    ProtocolCommSent = curComm;
#endif

    Bit x(false, PUBLIC);
    fmt::print("finish.\n", x.reveal());
}

void hist_min(
    std::vector<ShareA> &sharing_hist_vec,
    int bw_hist,
    int len_hist,
    uint64_t data_space_left, // 明文空间左边界
    uint64_t data_space_step, // 明文空间增量
    int bw_numeric            // 返回值的bit长度，也是明文空间的bit长度
)
{

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


// test 02
#ifdef LOG_LAYERWISE
    INIT_TIMER;
    INIT_ALL_IO_DATA_SENT;
#endif

    uint64_t mask_hist = my::mask_gen_uint64_t(bw_hist);
    auto sharing_v_vec = std::vector<ShareA>(len_hist);
    sharing_v_vec[0] = sharing_hist_vec[0];
    for (int i = 1; i < len_hist; i++)
    {
        sharing_v_vec[i] = (sharing_v_vec[i - 1] + sharing_hist_vec[i]) & mask_hist;
    }
    std::cout << std::endl;
    setup_semi_honest(iopackArr[0]->io_GC, party, 1024 * size_1024);

    auto tmp = iopackArr[0]->get_comm();
    auto tmp2 = iopackArr[0]->get_rounds();

    auto shareA = std::make_unique<Integer[]>(len_hist);
    auto shareB = std::make_unique<Integer[]>(len_hist);
    auto compare_res_bool = std::make_unique<Bit[]>(len_hist);
    // auto compare_res_Integer = std::make_unique<Integer[]>(len_hist);
    auto compare_res_share_A_Integer = std::make_unique<Integer[]>(len_hist);
    auto compare_res_share_B_Integer = std::make_unique<Integer[]>(len_hist);

    sci::PRG128 prg;
    // auto share_res_A = std::make_unique<uint64_t[]>(len_hist);
    // auto share_res_B = std::make_unique<uint64_t[]>(len_hist);
    // prg.random_data(share_res_A.get(), len_hist * sizeof(uint64_t));
    // for (int i = 0; i < len_hist; i++)
    // {
    //     share_res_A[i] = share_res_A[i] & mask_hist;
    // }

    auto zero = Integer(bw_hist, 0, PUBLIC);
    for (int i = 0; i < len_hist; i++)
    {
        shareA[i] = Integer(bw_hist, sharing_hist_vec[i], ALICE);
    }
    for (int i = 0; i < len_hist; i++)
    {
        shareB[i] = Integer(bw_hist, sharing_hist_vec[i], BOB);
    }
    for (int i = 0; i < len_hist; i++)
    {
        compare_res_bool[i] = ((shareA[i] + shareB[i]) == zero);
    }

    for (int i = 0; i < len_hist; i++)
    {
        compare_res_share_B_Integer[i] = If(
            compare_res_bool[i],
            Integer(bw_numeric, 1 , PUBLIC),
            Integer(bw_numeric, 0, PUBLIC));
    }

    Integer sharing_count_equal(bw_numeric, 0, PUBLIC);

    for (int i = 0; i < len_hist; i++)
    {
        sharing_count_equal = sharing_count_equal + compare_res_share_B_Integer[i];
    }

    uint64_t share_res_A;
    prg.random_data(&share_res_A, 1 * sizeof(uint64_t));
    
    uint64_t share_res_B = sharing_count_equal.reveal<uint64_t>(BOB);

    fmt::print(
        fg(fmt::color::red), "Equal->Share len_hist={}, GC: {} MiB, {} rounds\n",
        len_hist,
        (iopackArr[0]->get_comm() - tmp) / (1.0 * (1ULL << 20)),
        iopackArr[0]->get_rounds() - tmp2);
    fmt::print(fg(fmt::color::red), "AND Gates = {}\n", circ_exec->num_and());

    logger->info("Equal->Share len_hist={}, GC: {} MiB, {} rounds",
               len_hist,
               (iopackArr[0]->get_comm() - tmp) / (1.0 * (1ULL << 20)),
               iopackArr[0]->get_rounds() - tmp2);
    logger->info("AND Gates = {}", circ_exec->num_and());

// test 03
#ifdef LOG_LAYERWISE
    auto temp = TIMER_TILL_NOW;
    ProtocolTimeInMilliSec = temp;
    fmt::print(
        fg(fmt::color::red), "Time in sec for current Protocol = {}, ProtocolCommSent={}\n",
        (ProtocolTimeInMilliSec / 1000.0),
        ((ProtocolCommSent) / (1.0 * (1ULL << 20))));
    logger->info("Time in sec for current Protocol = {}, ProtocolCommSent={}",
            (ProtocolTimeInMilliSec / 1000.0),
            ((ProtocolCommSent) / (1.0 * (1ULL << 20))));
    uint64_t curComm;
    FIND_ALL_IO_TILL_NOW(curComm);
    ProtocolCommSent = curComm;
#endif

    Bit x(false, PUBLIC);
    fmt::print("finish.\n", x.reveal());
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
    amap.arg("sgc", size_1024, "gc batch size");
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

    log_filename = fmt::format("logs/my_gc_min_{}.txt", party);

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
               "-------------- {}-> min ----------------------\nnumber_of_data:{}, bw_numeric:{}, len_hist:{}, bw_hist:{}, left:{}, generate_min:{}, generate_max:{}, sgc={}\n",
               party_name, number_of_data, bw_numeric, len_hist, bw_hist, left, generate_min, generate_max,size_1024);
    logger->info("-------------- {}-> min ----------------------\nnumber_of_data:{}, bw_numeric:{}, len_hist:{}, bw_hist:{}, left:{}, generate_min:{}, generate_max:{}, sgc={}\n",
                 party_name, number_of_data, bw_numeric, len_hist, bw_hist, left, generate_min, generate_max,size_1024);

    StartComputation();

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

    hist_min(sharing_hist_vec, bw_hist, len_hist, left, step, bw_numeric);

    EndComputation();

    // test 06
    fmt::print(fg(fmt::color::red), "ProtocolTimeInMilliSec={}\n,ProtocolCommSent={} MiB\n", ProtocolTimeInMilliSec, ((ProtocolCommSent) / (1.0 * (1ULL << 20))));
    if (party == SERVER)
    {
        uint64_t ProtocolCommSentClient = 0;
        io->recv_data(&ProtocolCommSentClient, sizeof(uint64_t));
        fmt::print(fg(fmt::color::red), "ProtocolComm(sent+received) = {} MiB\n", ((ProtocolCommSent + ProtocolCommSentClient) / (1.0 * (1ULL << 20))));
    }
    else if (party == CLIENT)
    {
        io->send_data(&ProtocolCommSent, sizeof(uint64_t));
    }
}
