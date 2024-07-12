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
int32_t bitlength = 32;
spdlog::filename_t log_filename;

// test 01
uint64_t ProtocolTimeInMilliSec = 0;
uint64_t ProtocolCommSent = 0;

int size_1024 = 16;

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

void chi2(
    std::vector<ShareA> &sharing_hist_vec,
    int M, int L,
    int bw_hist,
    int K = -1)
{

    std::shared_ptr<spdlog::logger> logger;
    try
    {
        logger = spdlog::basic_logger_st<spdlog::async_factory>("chi2", log_filename);
    }
    catch (const spdlog::spdlog_ex &ex)
    {
        std::cout << "Log init failed: " << ex.what() << std::endl;
    }
    string party_name = (party == ALICE ? "ALICE" : "BOB");
    auto rounds_begin = iopackArr[0]->get_rounds();
    fmt::print(fg(fmt::color::green), "{}->rounds_begin= {}, bw_hist= {}, M= {}, L= {}\n", party_name, rounds_begin,bw_hist,M,L);
    logger->info("{}->rounds_begin= {}, bw_hist= {}, M= {}, L= {}", party_name, rounds_begin,bw_hist,M,L);


#ifdef LOG_LAYERWISE
    INIT_TIMER;
    INIT_ALL_IO_DATA_SENT;
#endif
    size_t len_hist = M * L;
    uint64_t mask_frequency = (bitlength == 64 ? -1 : ((1ULL << bitlength) - 1));

    std::vector<ShareA> sharing_R_mat(len_hist);
    std::vector<ShareA> sharing_P_mat(len_hist);

    generate_contingency_table(sharing_hist_vec, M, L, sharing_R_mat, sharing_P_mat);

    std::vector<ShareA> sharing_double_hist_vec = my::append_vectors(sharing_hist_vec, sharing_hist_vec);
    std::vector<ShareA> sharing_RP_vec = my::append_vectors(sharing_R_mat, sharing_P_mat);
    std::vector<ShareA> sharing_uv_vec(len_hist * 2);

    sci::PRG128 prg;

    setup_semi_honest(iopackArr[0]->io_GC, party, 1024 * size_1024);
    auto tmp = iopackArr[0]->get_comm();
    auto tmp2 = iopackArr[0]->get_rounds();

    Integer *sharing_double_hist_vec_0 = new Integer[len_hist * 2];
    Integer *sharing_double_hist_vec_1 = new Integer[len_hist * 2];
    Integer *sharing_RP_vec_0 = new Integer[len_hist * 2];
    Integer *sharing_RP_vec_1 = new Integer[len_hist * 2];
    Integer *UV_vec_Integer = new Integer[len_hist * 2];
    Integer *W_vec_Integer = new Integer[len_hist];

    block128 *compare_res_block = new block128[len_hist * 2];
    bool *compare_res_reveal = new bool[len_hist * 2];

    int W_Integer_bw = bw_hist * 3 + ceil(log2(len_hist));
    Integer W_Integer(bw_hist * 2, 0, PUBLIC);
    Integer K_Integer(W_Integer_bw, 0, PUBLIC);
    Integer shifted_1(W_Integer_bw, 1<<bw_hist, PUBLIC);

    ShareA sharing_chi2_A;
    prg.random_data(&sharing_chi2_A, sizeof(uint64_t));
    auto sharing_chi2_A_Integer = Integer(W_Integer_bw, sharing_chi2_A, ALICE);

    for (int i = 0; i < len_hist * 2; i++)
    {
        sharing_double_hist_vec_0[i] = Integer(bw_hist * 2, sharing_double_hist_vec[i], ALICE);
        sharing_RP_vec_0[i] = Integer(bw_hist * 2, sharing_RP_vec[i], ALICE);
    }

    for (int i = 0; i < len_hist * 2; i++)
    {
        sharing_double_hist_vec_1[i] = Integer(bw_hist * 2, sharing_double_hist_vec[i], BOB);
        sharing_RP_vec_1[i] = Integer(bw_hist * 2, sharing_RP_vec[i], BOB);
    }
    // Integer FourA = Integer(bw_hist, 4, ALICE);
    // Integer FourB = Integer(bw_hist, 4, BOB);
    for (int i = 0; i < len_hist * 2; i++)
    {
        UV_vec_Integer[i] = ((sharing_double_hist_vec_0[i] + sharing_double_hist_vec_1[i]) << bw_hist) / (sharing_RP_vec_0[i] + sharing_RP_vec_1[i]);
    }

    for (int i = 0; i < len_hist; i++)
    {
        W_vec_Integer[i] = UV_vec_Integer[i] * UV_vec_Integer[i + len_hist];
    }
    for (int i = 0; i < len_hist; i++)
    {
        W_Integer = W_Integer + W_vec_Integer[i];
    }

    Integer chi2_Integer = (W_Integer.resize(W_Integer_bw)- shifted_1)*K_Integer;
    
    // fmt::print(fg(fmt::color::red), "Chi2 len_hist={}, GC: {} MiB, {} rounds\n", len_hist,
    //            (iopackArr[0]->get_comm() - tmp) / (1.0 * (1ULL << 20)), iopackArr[0]->get_rounds() - tmp2);

    Integer shareing_chi2_B_Integer = chi2_Integer- sharing_chi2_A_Integer;
    ShareA sharing_chi2_B = shareing_chi2_B_Integer.reveal<uint64_t>(BOB);

// test 03
#ifdef LOG_LAYERWISE
    auto temp = TIMER_TILL_NOW;
    ProtocolTimeInMilliSec = temp;
    uint64_t curComm;
    FIND_ALL_IO_TILL_NOW(curComm);
    ProtocolCommSent = curComm;

    fmt::print(fg(fmt::color::red), "ProtocolTimeInMilliSec= {} ms, ProtocolCommSent= {} MiB\n", ProtocolTimeInMilliSec, ((ProtocolCommSent) / (1.0 * (1ULL << 20))));

    fmt::print(fg(fmt::color::red), "reveal len_hist={}, GC: {} MiB, {} rounds\n", len_hist,
               (iopackArr[0]->get_comm() - tmp) / (1.0 * (1ULL << 20)), iopackArr[0]->get_rounds() - tmp2);
               

    logger->info("ProtocolTimeInMilliSec= {} ms, ProtocolCommSent= {} MiB", ProtocolTimeInMilliSec, ((ProtocolCommSent) / (1.0 * (1ULL << 20))));
    logger->info("reveal len_hist={}, GC: {} MiB, {} rounds", len_hist,
               (iopackArr[0]->get_comm() - tmp) / (1.0 * (1ULL << 20)), iopackArr[0]->get_rounds() - tmp2);
    // if (party == SERVER)
    // {
    //     uint64_t ProtocolCommSentClient = 0;
    //     io->recv_data(&ProtocolCommSentClient, sizeof(uint64_t));
    //     fmt::print(fg(fmt::color::red), "ProtocolComm(sent+received) = {} MiB\n", ((ProtocolCommSent + ProtocolCommSentClient) / (1.0 * (1ULL << 20))));
    // }
    // else if (party == CLIENT)
    // {
    //     io->send_data(&ProtocolCommSent, sizeof(uint64_t));
    // }

#endif

    Bit x(false, PUBLIC);
    fmt::print("finish.\n", x.reveal());

    delete[] sharing_double_hist_vec_0;
    delete[] sharing_double_hist_vec_1;
    delete[] sharing_RP_vec_0;
    delete[] sharing_RP_vec_1;
    delete[] UV_vec_Integer;
    delete[] compare_res_block;
    delete[] compare_res_reveal;

}

int main(int argc, char **argv)
{
    int left = -1;
    int step = -1;
    int number_of_data = -1;
    uint64_t generate_min = -1;
    uint64_t generate_max = -1;
    int bw_hist = -1;
    int bw_numeric = -1;
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
    amap.arg("bh", bw_hist, "Bitwidth of histgram");
    amap.arg("min", generate_min, "生成的数据的最小值");
    amap.arg("max", generate_max, "生成的数据的最大值");
    amap.arg("nd", number_of_data, "生成的数据的个数");
    amap.arg("M", M, "M");
    amap.arg("L", L, "L");
    amap.arg("bn", bw_numeric, "Bitwidth of numeric data, 在chi2中应当是ceil(log2(number_of_data))的2倍以上");
    amap.arg("sgc", size_1024, "gc batch size");
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
    my::set_value_int(bw_hist, ceil(log2(number_of_data)) + 1);
    my::set_value_int(bw_numeric, ceil(log2(left + (len_hist - 1) * step)));

    log_filename = fmt::format("logs/my_gc_chi2_{}.txt", party);

    dim_frequency_global = len_hist;
    num_data_global = number_of_data;

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
               "-------------- {}-> chi2 ----------------------\nnumber_of_data:{}, bw_numeric:{}, len_hist:{}, bw_hist:{},  left:{}, generate_min:{}, generate_max:{}, sgc={}\n",
               party_name, number_of_data, bw_numeric, len_hist, bw_hist, left, generate_min, generate_max,size_1024);
    logger->info("-------------- {}-> chi2 ----------------------\nnumber_of_data:{}, bw_numeric:{}, len_hist:{}, bw_hist:{},  left:{}, generate_min:{}, generate_max:{}, sgc={}\n",
                 party_name, number_of_data, bw_numeric, len_hist, bw_hist, left, generate_min, generate_max,size_1024);

    StartComputation();
    std::vector<uint64_t> hist_vec_actual_ALICE;
    std::vector<uint64_t> sharing_hist_vec;

    auto iopackMain = new IOPack(party, port - 1, address);

    if (party == ALICE)
    {
        auto numeric_data_vec = my::generate_numeric_data_vector(number_of_data, bw_numeric, generate_min, generate_max);
        hist_vec_actual_ALICE = my::numeric_data_to_hist(len_hist, left, step, numeric_data_vec, number_of_data);
        sharing_hist_vec = my::generate_and_send_sharing_vector_plaintext(iopackMain, hist_vec_actual_ALICE, bw_hist);
    }
    else if (party == BOB)
    {
        sharing_hist_vec = my::recv_sharing_vector_plaintext_vector(iopackMain, len_hist);
    }
    delete iopackMain;

    chi2(sharing_hist_vec, M, L, bw_hist);
    EndComputation();

}
