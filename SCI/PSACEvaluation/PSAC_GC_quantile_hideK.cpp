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
// string address = "172.29.20.64";
string address = "127.0.0.1";
spdlog::filename_t log_filename;
int32_t bitlength = 32;

// test 01
uint64_t ProtocolTimeInMilliSec = 0;
uint64_t ProtocolCommSent = 0;

int size_1024 = 16;

struct QuantileInfo
{
    uint64_t p_prime;
    int l_prime;
    double percentile;
    int num_of_data;
    int index;
    double gamma;
};

struct QuantileArithmetic
{
    uint64_t left;
    uint64_t right;
} quantile_actual;
QuantileInfo quantile_number(uint64_t p_prime, int l_prime, int num_of_data)
{
    QuantileInfo qinfo{};
    qinfo.p_prime = p_prime;
    qinfo.l_prime = l_prime;
    qinfo.percentile = double(p_prime) / pow(2, l_prime);
    qinfo.num_of_data = num_of_data;
    qinfo.index = floor((num_of_data - 1) * qinfo.percentile) + 1;
    qinfo.gamma = (num_of_data - 1) * qinfo.percentile - floor((num_of_data - 1) * qinfo.percentile);
    return qinfo;
}

ShareA hist_quantile(
    uint64_t p_prime, int l_prime,
    std::vector<ShareA> &sharing_hist_vec,
    uint64_t data_space_left,
    uint64_t data_space_step,
    int bw_hist)
{

    int bw_compare = bw_hist - l_prime;
    int len_hist = sharing_hist_vec.size();

    std::shared_ptr<spdlog::logger> logger;
    try
    {
        logger = spdlog::basic_logger_st<spdlog::async_factory>(fmt::format("quantile_{}/{}", p_prime, int(pow(2, l_prime))), log_filename);
    }
    catch (const spdlog::spdlog_ex &ex)
    {
        std::cout << "Log init failed: " << ex.what() << std::endl;
    }
    string party_name = (party == ALICE ? "ALICE" : "BOB");
    auto rounds_begin = iopackArr[0]->get_rounds();

    fmt::print(fg(fmt::color::green), "{}->rounds_begin= {}, bw_hist= {}, bw_compare= {}\n", party_name, rounds_begin, bw_hist, bw_compare);
    logger->info("{}->rounds_begin= {}, bw_hist= {}, bw_compare= {}", party_name, rounds_begin, bw_hist, bw_compare);

// test 02
#ifdef LOG_LAYERWISE
    INIT_TIMER;
    INIT_ALL_IO_DATA_SENT;
#endif
    uint64_t mask_hist = my::mask_gen_uint64_t(bw_hist);
    uint64_t mask_compare = my::mask_gen_uint64_t(bw_compare);

    sci::PRG128 prg;

    setup_semi_honest(iopackArr[0]->io_GC, party, 1024 * size_1024);
    auto tmp = iopackArr[0]->get_comm();
    auto tmp2 = iopackArr[0]->get_rounds();

    auto share_1 = my::shareA_const(party, 1);
    ShareA sharinf_K_prime = std::reduce(std::execution::par_unseq, sharing_hist_vec.begin(), sharing_hist_vec.end(), 0) & mask_hist;
    ShareA j_tmp = (sharinf_K_prime - share_1) * p_prime & mask_hist;
    ShareA j_tmp_out;

    auto sharing_v_vec = std::vector<ShareA>(len_hist);
    sharing_v_vec[0] = sharing_hist_vec[0];
    for (int i = 1; i < len_hist; i++)
    {
        sharing_v_vec[i] = (sharing_v_vec[i - 1] + sharing_hist_vec[i]) & mask_compare;
    }

    ShareA sharing_gamma = (j_tmp - (j_tmp_out << l_prime) & mask_hist) & mask_hist;

    Integer Integer_K_share_B = Integer(bw_hist, sharinf_K_prime, BOB);
    Integer Integer_gamma_B = Integer(bw_hist, sharing_gamma, BOB);

    Integer Integer_K_share_A = Integer(bw_hist, sharinf_K_prime, ALICE);
    Integer Integer_gamma_A = Integer(bw_hist, sharing_gamma, ALICE);

    Integer Integer_gamma = Integer_gamma_A + Integer_gamma_B;
    Integer Integer_j = ((Integer_K_share_A + Integer_K_share_B) >> l_prime) + Integer(bw_compare, 1, PUBLIC);

    auto shareA = std::make_unique<Integer[]>(len_hist);
    auto shareB = std::make_unique<Integer[]>(len_hist);
    Integer Integer_Qp_share_A;
    uint64_t Qp_share_A;

    for (int i = 0; i < len_hist; i++)
    {
        shareA[i] = Integer(bw_hist, sharing_hist_vec[i], ALICE);
    }
    for (int i = 0; i < len_hist; i++)
    {
        shareB[i] = Integer(bw_hist, sharing_hist_vec[i], BOB);
    }
    auto compare_res_Bit = std::make_unique<Bit[]>(len_hist);
    auto equal_res_Bit = std::make_unique<Bit[]>(len_hist);
    for (int i = 0; i < len_hist; i++)
    {
        compare_res_Bit[i] = ((shareA[i] + shareB[i]) < Integer_j);
    }
    for (int i = 0; i < len_hist; i++)
    {
        equal_res_Bit[i] = ((shareA[i] + shareB[i]) == Integer_j);
    }

    fmt::print(fg(fmt::color::red), "Compare->Bit len_hist={}, GC: {} MiB, {} rounds\n", len_hist,
               (iopackArr[0]->get_comm() - tmp) / (1.0 * (1ULL << 20)), iopackArr[0]->get_rounds() - tmp2);

    // auto share_compare_res_A = std::make_unique<uint64_t[]>(len_hist);
    // auto share_equal_res_A = std::make_unique<uint64_t[]>(len_hist);
    // prg.random_data(share_compare_res_A.get(), len_hist * sizeof(uint64_t));
    // prg.random_data(share_equal_res_A.get(), len_hist * sizeof(uint64_t));

    auto compare_res_Integer = std::make_unique<Integer[]>(len_hist);
    auto equal_res_Integer = std::make_unique<Integer[]>(len_hist);
    for (int i = 0; i < len_hist; i++)
    {
        compare_res_Integer[i] = If(
            compare_res_Bit[i],
            Integer(bw_hist, 1, PUBLIC),
            Integer(bw_hist, 0, PUBLIC));
        equal_res_Integer[i] = If(
            equal_res_Bit[i],
            Integer(bw_hist, 1, PUBLIC),
            Integer(bw_hist, 0, PUBLIC));
    }
    // for (int i = 0; i < len_hist; i++)
    // {
    // }

    fmt::print(fg(fmt::color::red), "Compare->Integer len_hist={}, GC: {} MiB, {} rounds\n", len_hist,
               (iopackArr[0]->get_comm() - tmp) / (1.0 * (1ULL << 20)), iopackArr[0]->get_rounds() - tmp2);

    // auto share_compare_res_B = std::make_unique<uint64_t[]>(len_hist);
    // auto share_equal_res_B = std::make_unique<uint64_t[]>(len_hist);

    Integer count_less = Integer(bw_hist, 0, PUBLIC);
    Integer count_eq = Integer(bw_hist, 0, PUBLIC);

    for (int i = 0; i < len_hist; i++)
    {
        // share_compare_res_B[i] = compare_res_Integer[i].reveal<uint64_t>(BOB);
        count_less = count_less + compare_res_Integer[i];
        count_eq = count_eq + equal_res_Integer[i];
    }
    Integer W = Integer_gamma * count_eq;
    Integer Q = ((Integer(bw_hist, data_space_left, PUBLIC) + count_less) << l_prime) + W;

    uint64_t share_res_A;
    prg.random_data(&share_res_A, 1 * sizeof(uint64_t));
    uint64_t share_res_B = Q.reveal<uint64_t>(BOB);

    fmt::print(
        fg(fmt::color::red),
        "Compare->Share len_hist={}, GC: {} MiB, {} rounds\n",
        len_hist,
        (iopackArr[0]->get_comm() - tmp) / (1.0 * (1ULL << 20)),
        iopackArr[0]->get_rounds() - tmp2);
    fmt::print(fg(fmt::color::red), "AND Gates = {}\n", circ_exec->num_and());

    logger->info("Compare->Share len_hist={}, GC: {} MiB, {} rounds",
                 len_hist,
                 (iopackArr[0]->get_comm() - tmp) / (1.0 * (1ULL << 20)),
                 iopackArr[0]->get_rounds() - tmp2);
    logger->info("AND Gates = {}", circ_exec->num_and());

// test 03
#ifdef LOG_LAYERWISE
    auto temp = TIMER_TILL_NOW;
    ProtocolTimeInMilliSec = temp;
    fmt::print(fg(fmt::color::red), "Time in sec for current Protocol = {}, ProtocolCommSent={}\n", (ProtocolTimeInMilliSec / 1000.0), ((ProtocolCommSent) / (1.0 * (1ULL << 20))));
    
    logger->info("Time in sec for current Protocol = {}, ProtocolCommSent={}",
            (ProtocolTimeInMilliSec / 1000.0),
            ((ProtocolCommSent) / (1.0 * (1ULL << 20))));
    uint64_t curComm;
    FIND_ALL_IO_TILL_NOW(curComm);
    ProtocolCommSent = curComm;
#endif

    Bit x(false, PUBLIC);
    fmt::print("finish.\n", x.reveal());

    if (party == ALICE)
    {
        return share_res_A;
    }
    else
    {
        return share_res_B;
    }
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
    // double percentile = -1.0;
    uint64_t p_prime = -1;
    int l_prime = -1;

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
    // amap.arg("pe", percentile, "percentile");
    amap.arg("pr", p_prime, "p_prime");
    amap.arg("lr", l_prime, "l_prime");
    amap.arg("sgc", size_1024, "gc batch size");
    amap.parse(argc, argv);

    my::set_value_int(len_hist, 10000);
    my::set_value_int(left, 100);
    my::set_value_int(step, 1);

    my::set_value_int(number_of_data, 4000);
    my::set_value_uint64_t(generate_min, left + len_hist * step * 0.02);
    my::set_value_uint64_t(generate_max, left + (len_hist - 1) * step - len_hist * step * 0.02);

    // my::set_value_double(percentile, 0.1);
    my::set_value_uint64_t(p_prime, 1);
    my::set_value_int(l_prime, 2);
    int bw_hist1 = ceil(log2(number_of_data)) + l_prime + 1;
    int bw_hist2 = ceil(log2(left + (len_hist - 1) * step)) + l_prime;
    my::set_value_int(bw_hist, bw_hist1 > bw_hist2 ? bw_hist1 : bw_hist2);
    my::set_value_int(bw_numeric, ceil(log2(left + (len_hist - 1) * step)));

    log_filename = fmt::format("logs/my_gc_quantile_hideK_{}.txt", party);

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
    auto qinfo = quantile_number(p_prime, l_prime, number_of_data);
    fmt::print(fg(fmt::color::blue),
               "-------------- {}-> quantile_hideK ----------------------\nnumber_of_data:{}, bw_numeric:{}, len_hist:{}, bw_hist:{}, left:{}, generate_min:{}, generate_max:{}, percentile:{}, p_prime= {}, l_prime= {}, bw_hist1= {},bw_hist2= {}, sgc={}\n",
               party_name, number_of_data, bw_numeric, len_hist, bw_hist, left, generate_min, generate_max, qinfo.percentile, p_prime, l_prime, bw_hist1, bw_hist2,size_1024);
    logger->info(
        "-------------- {}-> quantile_hideK ----------------------\nnumber_of_data:{}, bw_numeric:{}, len_hist:{}, bw_hist:{}, left:{}, generate_min:{}, generate_max:{}, percentile:{}, p_prime= {}, l_prime= {}, bw_hist1= {},bw_hist2= {}, sgc={}\n",
        party_name, number_of_data, bw_numeric, len_hist, bw_hist, left, generate_min, generate_max, qinfo.percentile, p_prime, l_prime, bw_hist1, bw_hist2,size_1024);

    StartComputation();
    std::vector<uint64_t> sharing_hist_vec;

    auto iopackMain = new IOPack(party, port - 1, address);

    if (party == ALICE)
    {
        auto numeric_data_vec = my::generate_numeric_data_vector(number_of_data, bw_numeric, generate_min, generate_max);
        auto hist_vec = my::numeric_data_to_hist(len_hist, left, step, numeric_data_vec, number_of_data);
        std::vector<uint64_t> data_accumulate_vec(len_hist);
        data_accumulate_vec[0] = hist_vec[0];
        for (int i = 1; i < len_hist; i++)
        {
            data_accumulate_vec[i] = data_accumulate_vec[i - 1] + hist_vec[i];
        }
        sharing_hist_vec = my::generate_and_send_sharing_vector_plaintext(iopackArr[0], hist_vec, bw_hist);

#ifdef _VERIFY_
        vector<uint64_t> data_vec(numeric_data_vec);
        std::sort(data_vec.begin(), data_vec.end());
        quantile_actual.left = data_vec[qinfo.index - 1];
        quantile_actual.right = data_vec[qinfo.index];
        fmt::print(fg(fmt::color::blue), "Line:{}, quantile_actual_left={}, quantile_actual_right={}\n", __LINE__, quantile_actual.left, quantile_actual.right);
#endif
    }
    else if (party == BOB)
    {
        sharing_hist_vec = my::recv_sharing_vector_plaintext_vector(iopackArr[0], len_hist);
    }

    ShareA sharing_Q_p = hist_quantile(p_prime, l_prime, sharing_hist_vec, left, step, bw_hist);

    EndComputation();

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

    /************** Verification ****************/
#ifdef _VERIFY_
    std::vector<ShareA> sharing_Q_p_vec(1, sharing_Q_p);
    auto Q_p_vec = my::reconstruct_uint64_t(party, sci::ALICE, iopackArr, sharing_Q_p_vec, bw_hist);
    if (party == ALICE)
    {
        double Q_p_double = Q_p_vec[0] / pow(2, l_prime);
        double Q_p_double_actual = (1 - qinfo.gamma) * quantile_actual.left + qinfo.gamma * quantile_actual.right;
        fmt::print(fg(fmt::color::blue), "Line:{}, Q_p= {}, Q_p_double= {}, Q_p_double_actual= {}\n", __LINE__, Q_p_vec, Q_p_double, Q_p_double_actual);
        assert(Q_p_double_actual == Q_p_double);
    }
#endif
}
