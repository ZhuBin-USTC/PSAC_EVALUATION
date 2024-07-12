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

// bw_hist_part1 = max[ceil(log2(number_of_data)), ceil(log2((len_hist - 1) * step)+1)];
// bw_hist_part2 = max[ceil(log2(k_double_prime)), l_double_prime]
// bw_hist = max[bw_hist_part1 + bw_hist_part2 + 2 + 1, ceil(log2(left + (len_hist - 1) * step)) + 2];

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
} quantile_1_4_actual, quantile_3_4_actual;
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
    uint64_t mask_hist = my::mask_gen_uint64_t(bw_hist);
    uint64_t mask_compare = my::mask_gen_uint64_t(bw_compare);

#ifdef _PERFORMANCE_
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

    uint64_t last_rounds = rounds_begin;

    uint64_t total_comm = 0;
    uint64_t thread_comm[num_threads];
    for (int i = 0; i < num_threads; i++)
    {
        thread_comm[i] = iopackArr[i]->get_comm();
    }
    auto start = sci::clock_start();
#endif

    auto share_1 = my::shareA_const(party, 1);
    ShareA sharinf_K_prime = std::reduce(std::execution::par_unseq, sharing_hist_vec.begin(), sharing_hist_vec.end(), 0) & mask_hist;
    ShareA tmp = (sharinf_K_prime - share_1) * p_prime & mask_hist;
    ShareA tmp_out;

    auto linearOT_ptr = std::make_unique<LinearOT>(party, iopackArr[0], otpackArr[0]);
    // 交互1
    // my::truncate_thread(TRUNCATE_RED_THEN_EXT, party, iopackArr, otpackArr, 0, 1, &tmp, &tmp_out, l_prime, bw_hist, false);

    // 16bit
    // linearOT_ptr->trunc->truncate_red_then_ext(1,&tmp, &tmp_out,l_prime,bw_hist,false); //10
    // linearOT_ptr->trunc->truncate(1,&tmp, &tmp_out,l_prime,bw_hist,false); //12
    // linearOT_ptr->trunc->div_pow2(1,&tmp, &tmp_out,l_prime,bw_hist,false); //12
    linearOT_ptr->trunc->truncate_and_reduce(1, &tmp, &tmp_out, l_prime, bw_hist); // 4，没有extern

#ifdef _PERFORMANCE_
    uint64_t truncate_rounds = iopackArr[0]->get_rounds() - last_rounds;
    last_rounds = iopackArr[0]->get_rounds();
    fmt::print(fg(fmt::color::white_smoke), "{}->truncate rounds= {}, num rounds= {}\n", party_name, truncate_rounds, last_rounds);
    logger->info("{}->truncate rounds= {}, num rounds= {}", party_name, truncate_rounds, last_rounds);
#endif

    ShareA j_A = (tmp_out + share_1) & mask_compare;
    ShareA gamma_A = (tmp - (tmp_out << l_prime) & mask_hist) & mask_hist;
    // gamma shifted left by l'.

#ifdef _CHECK_
    fmt::print(fg(fmt::color::red), "sharinf_K_prime= {}, tmp= {}, tmp_out= {}, j_A ={}, gamma_A= {}\n", sharinf_K_prime, tmp, tmp_out, j_A, gamma_A);
#endif

    auto sharing_v_vec = std::vector<ShareA>(len_hist);
    sharing_v_vec[0] = sharing_hist_vec[0];
    for (int i = 1; i < len_hist; i++)
    {
        sharing_v_vec[i] = (sharing_v_vec[i - 1] + sharing_hist_vec[i]) & mask_compare;
    }

    auto less_vec_B = std::vector<ShareB>(len_hist); // 比较大小返回值向量
    auto less_vec_A = std::vector<ShareA>(len_hist);
    auto eq_vec_B = std::vector<ShareB>(len_hist); // 比较大小返回值向量
    auto eq_vec_A = std::vector<ShareA>(len_hist);

    auto sharing_cmp_vec = std::vector<ShareA>(len_hist);
    std::transform(std::execution::par_unseq, sharing_v_vec.begin(), sharing_v_vec.end(), sharing_cmp_vec.begin(), [j_A](ShareA a)
                   { return a - j_A; });

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

        // 交互2
        compare_arithmetic_threads[tid] = std::thread(
            // my::compare_arithmetic_thread, GREATERwithEQUAL, party, iopackArr, otpackArr, tid,
            my::compare_arithmetic_thread, LESSwithEQUAL, party, iopackArr, otpackArr, tid,
            sharing_cmp_vec.data() + offset, 0, num_ops, bw_compare, bw_hist,
            less_vec_B.data() + offset, eq_vec_B.data() + offset, less_vec_A.data() + offset, eq_vec_A.data() + offset);
    }
    for (int tid = 0; tid < num_threads; ++tid)
    {
        compare_arithmetic_threads[tid].join();
    }

#ifdef _PERFORMANCE_
    uint64_t compare_rounds = iopackArr[0]->get_rounds() - last_rounds;
    last_rounds = iopackArr[0]->get_rounds();
    fmt::print(fg(fmt::color::white_smoke), "{}->compare_rounds= {}, num rounds= {}\n", party_name, compare_rounds, last_rounds);
    logger->info("{}->compare_rounds= {}, num rounds= {}", party_name, compare_rounds, last_rounds);
#endif

    uint64_t count_less_eq_A_res[2] = {0}; // 0是less，1是equal
    count_less_eq_A_res[0] = std::reduce(std::execution::par_unseq, less_vec_A.begin(), less_vec_A.begin() + len_hist, 0) & mask_hist;
    count_less_eq_A_res[1] = std::reduce(std::execution::par_unseq, eq_vec_A.begin(), eq_vec_A.begin() + len_hist, 0) & mask_hist;

    // 交互3
    ShareA sharing_w;
    linearOT_ptr->hadamard_product(1, &gamma_A, count_less_eq_A_res + 1, &sharing_w, bw_hist, bw_hist, bw_hist, false, false, MultMode::None, nullptr, nullptr);

    ShareA sharing_Q_p = (((my::shareA_const(party, data_space_left) + count_less_eq_A_res[0]) << l_prime) + sharing_w) & mask_hist;

#ifdef _VERIFY_
    fmt::print(fg(fmt::color::blue), "sharing_Q_p = {}\n", sharing_Q_p);
#endif

    // ShareA count_less_eq_A_num[2]; // 1的数量是max结果的index的sharing
    // linearOT_ptr->xt->z_extend(1, &sharing_Q_p, &sharing_Q_p, bw_hist, bw_numeric);

#ifdef _PERFORMANCE_
    long long t = time_from(start);
    for (int i = 0; i < num_threads; i++)
    {
        thread_comm[i] = iopackArr[i]->get_comm() - thread_comm[i];
        total_comm += thread_comm[i];
    }

    uint64_t hadamard_product_rounds = iopackArr[0]->get_rounds() - last_rounds;
    last_rounds = iopackArr[0]->get_rounds();
    fmt::print(fg(fmt::color::white_smoke), "{}->hadamard_product_rounds= {}, num rounds= {}\n", party_name, hadamard_product_rounds, last_rounds);
    logger->info("{}->hadamard_product_rounds= {}, num rounds= {}", party_name, hadamard_product_rounds, last_rounds);

    fmt::print(fg(fmt::color::white_smoke), "{}->Number of dim/s:\t{}\n", party_name, (double(len_hist) / t) * 1e6);
    fmt::print(fg(fmt::color::white_smoke), "{}->Computing Time\t{} ms\n", party_name, t / (1000.0));
    fmt::print(fg(fmt::color::white_smoke), "{}->Computing Bytes Sent\t{} bytes, {} MB\n", party_name, total_comm, total_comm / (1024.0 * 1024.0));
    logger->info("{}->Number of dim/s:\t{}", party_name, (double(len_hist) / t) * 1e6);
    logger->info("{}->Computing Time\t{} ms", party_name, t / (1000.0));
    logger->info("{}->Computing Bytes Sent\t{} bytes, {} MB", party_name, total_comm, total_comm / (1024.0 * 1024.0));
#endif
    spdlog::drop_all();
    return sharing_Q_p;
}

std::vector<ShareB> outlier_detect(
    std::vector<ShareA> &sharing_hist_vec,
    uint64_t data_space_left,
    uint64_t data_space_step,
    ShareA Q_3_4,
    ShareA Q_1_4,
    int bw_hist,
    // int bw_numeric, // Q_1_4, Q_3_4
    int quantile_shifted,
    uint64_t k_double_prime,
    int l_double_prime)
{
    auto mask_hist = my::mask_gen_uint64_t(bw_hist);
    // auto mask_numeric = my::mask_gen_uint64_t(bw_numeric);

    ShareA sharing_I = k_double_prime * (Q_3_4 - Q_1_4);
    ShareA sharing_x_H = (((Q_3_4 - (my::shareA_const(party, data_space_left) << quantile_shifted)) << l_double_prime) + sharing_I) & mask_hist;
    ShareA sharing_x_L = (((Q_1_4 - (my::shareA_const(party, data_space_left) << quantile_shifted)) << l_double_prime) - sharing_I) & mask_hist;

#ifdef _VERIFY_
    fmt::print(fg(fmt::color::blue), "LINE= {}, sharing_x_H = {}, sharing_x_L = {}\n", __LINE__, sharing_x_H, sharing_x_L);
#endif

#ifdef _PERFORMANCE_
    std::shared_ptr<spdlog::logger> logger;
    try
    {
        logger = spdlog::basic_logger_st<spdlog::async_factory>("outlier_detect", log_filename);
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

    auto len_hist = sharing_hist_vec.size();
    std::vector<ShareA> sharing_x_H_vec(len_hist, sharing_x_H);
    std::vector<ShareA> sharing_x_L_vec(len_hist, sharing_x_L);
    std::vector<ShareA> sharing_s(len_hist, 0);

    if (party == ALICE)
    {
        uint64_t start_shifted = 0;
        uint64_t step_shifted = data_space_step << (l_double_prime + quantile_shifted);
        std::generate(sharing_s.begin(), sharing_s.end(), [start_shifted, step_shifted]()
                      { static uint64_t i = start_shifted-step_shifted; return i+=step_shifted; });
        std::transform(std::execution::par_unseq,
                       sharing_s.begin(), sharing_s.end(), sharing_s.begin(), [mask_hist](ShareA e)
                       { return e & mask_hist; });
    }

    std::vector<ShareA> sharing_x_H_append_s(sharing_x_H_vec);
    std::vector<ShareA> sharing_s_append_x_L(sharing_s);

    sharing_x_H_append_s.insert(sharing_x_H_append_s.end(), sharing_s.begin(), sharing_s.end());
    sharing_s_append_x_L.insert(sharing_s_append_x_L.end(), sharing_x_L_vec.begin(), sharing_x_L_vec.end());

    std::vector<ShareA> sharing_cmp_vec(sharing_x_H_append_s);
    std::transform(std::execution::par_unseq,
                   sharing_cmp_vec.begin(), sharing_cmp_vec.end(), sharing_s_append_x_L.begin(),
                   sharing_cmp_vec.begin(), std::minus<ShareA>());
    std::transform(std::execution::par_unseq,
                   sharing_cmp_vec.begin(), sharing_cmp_vec.end(),
                   sharing_cmp_vec.begin(), [mask_hist](ShareA e)
                   { return e & mask_hist; });

    std::vector<ShareB> sharing_V_H_append_V_L(2 * len_hist);
    std::vector<ShareB> sharing_V(len_hist);

    std::thread compare_threads[num_threads];
    int chunk_size = 2 * len_hist / num_threads;
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
        // 交互2
        compare_threads[tid] = std::thread(
            my::compare_thread, LESS, party, iopackArr, otpackArr, tid,
            sharing_cmp_vec.data() + offset, 0, num_ops, bw_hist, sharing_V_H_append_V_L.data() + offset, nullptr);
    }
    for (int tid = 0; tid < num_threads; ++tid)
    {
        compare_threads[tid].join();
    }

    std::transform(std::execution::par_unseq,
                   sharing_V_H_append_V_L.begin(), sharing_V_H_append_V_L.begin() + len_hist,
                   sharing_V_H_append_V_L.begin() + len_hist,
                   sharing_V.begin(), [](auto e1, auto e2)
                   { return (e1 + e2) % 2; });

#ifdef _PERFORMANCE_
    long long t = time_from(start);
    for (int i = 0; i < num_threads; i++)
    {
        thread_comm[i] = iopackArr[i]->get_comm() - thread_comm[i];
        total_comm += thread_comm[i];
    }

    uint64_t compare_rounds = iopackArr[0]->get_rounds() - last_rounds;
    last_rounds = iopackArr[0]->get_rounds();
    fmt::print(fg(fmt::color::white_smoke), "{}->compare_rounds= {}, num rounds= {}\n", party_name, compare_rounds, last_rounds);
    logger->info("{}->compare_rounds= {}, num rounds= {}", party_name, compare_rounds, last_rounds);

    /**** Process & Write Benchmarking Data *****/
    fmt::print(fg(fmt::color::white_smoke), "{}->Number of dim/s:\t{}\n", party_name, (double(len_hist) / t) * 1e6);
    fmt::print(fg(fmt::color::white_smoke), "{}->Computing Time\t{} ms\n", party_name, t / (1000.0));
    fmt::print(fg(fmt::color::white_smoke), "{}->Computing Bytes Sent\t{} bytes, {} MB\n", party_name, total_comm, total_comm / (1024.0 * 1024.0));
    logger->info("{}->Number of dim/s:\t{}", party_name, (double(len_hist) / t) * 1e6);
    logger->info("{}->Computing Time\t{} ms", party_name, t / (1000.0));
    logger->info("{}->Computing Bytes Sent\t{} bytes, {} MB", party_name, total_comm, total_comm / (1024.0 * 1024.0));
#endif

#ifdef _CHECK_
    auto x_H_append_s = my::reconstruct_uint64_t(party, sci::ALICE, iopackArr, sharing_x_H_append_s, bw_hist);
    auto s_append_x_L = my::reconstruct_uint64_t(party, sci::ALICE, iopackArr, sharing_s_append_x_L, bw_hist);
    auto cmp_vec = my::reconstruct_uint64_t(party, sci::ALICE, iopackArr, sharing_cmp_vec, bw_hist);
    auto V_H_append_V_L = my::reconstruct_uint8_t(party, sci::ALICE, iopackArr, sharing_V_H_append_V_L);
    fmt::print(fg(fmt::color::purple), "{}", x_H_append_s);
    fmt::print(fg(fmt::color::purple), "Line= {}\n", __LINE__);
    if (party == ALICE)
    {
        for (int i = 0; i < len_hist * 2; i++)
        {
            if (i % 10 == 0)
            {
                fmt::print(fg(fmt::color::purple), "i={}, x_H_append_s={}, s_append_x_L={}, cmp_vec={}, V_H_append_V_L={}\n", i, x_H_append_s[i], s_append_x_L[i], cmp_vec[i], V_H_append_V_L[i]);
            }
        }
    }
#endif
    /******************* Cleanup ****************/
    spdlog::drop_all();
    return sharing_V;
}

std::vector<ShareA> outlier_estimation(
    std::vector<ShareA> &sharing_hist_vec,
    std::vector<ShareB> &sharing_V,
    int bw_hist)
{
#ifdef _PERFORMANCE_
    std::shared_ptr<spdlog::logger> logger;
    try
    {
        logger = spdlog::basic_logger_st<spdlog::async_factory>("outlier_estimation", log_filename);
    }
    catch (const spdlog::spdlog_ex &ex)
    {
        std::cout << "Log init failed: " << ex.what() << std::endl;
    }
    string party_name = (party == ALICE ? "ALICE" : "BOB");
    auto rounds_begin = iopackArr[0]->get_rounds();
    fmt::print(fg(fmt::color::green), "{}->rounds_begin= {}\n", party_name, rounds_begin);
    logger->info("{}->rounds_begin= {}", party_name, rounds_begin);
    uint64_t last_rounds = rounds_begin;

    uint64_t total_comm = 0;
    uint64_t thread_comm[num_threads];
    for (int i = 0; i < num_threads; i++)
    {
        thread_comm[i] = iopackArr[i]->get_comm();
    }
    auto start = clock_start();
#endif

    auto num_rounds = iopackArr[0]->get_rounds();
    fmt::print(fg(fmt::color::green), "---------------------------\n outlier_estimation_begin num_rounds= {}\n", num_rounds);
    auto len_hist = sharing_hist_vec.size();

    std::vector<ShareA> sharing_hist_vec_star(len_hist);

    std::thread compare_threads[num_threads];
    int chunk_size = len_hist / num_threads;

    std::vector<ShareB> sharing_V_star(sharing_V);
    if (party == sci::ALICE)
    {
        std::transform(std::execution::par_unseq, sharing_V_star.begin(), sharing_V_star.end(), sharing_V_star.begin(), [](auto e)
                       { return (e ^ 1); });
    }

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
        // 交互2
        std::vector<ShareA> y;
        compare_threads[tid] = std::thread(
            my::mux_thread, party, iopackArr, otpackArr, tid, sharing_V_star.data() + offset, sharing_hist_vec.data() + offset, sharing_hist_vec_star.data() + offset, chunk_size, bw_hist, bw_hist);
    }
    for (int tid = 0; tid < num_threads; ++tid)
    {
        compare_threads[tid].join();
    }

#ifdef _PERFORMANCE_
    long long t = time_from(start);
    for (int i = 0; i < num_threads; i++)
    {
        thread_comm[i] = iopackArr[i]->get_comm() - thread_comm[i];
        total_comm += thread_comm[i];
    }

    uint64_t mux_rounds = iopackArr[0]->get_rounds() - last_rounds;
    last_rounds = iopackArr[0]->get_rounds();
    fmt::print(fg(fmt::color::white_smoke), "{}->mux_rounds= {}, num rounds= {}\n", party_name, mux_rounds, last_rounds);
    logger->info("{}->mux_rounds= {}, num rounds= {}", party_name, mux_rounds, last_rounds);

    /**** Process & Write Benchmarking Data *****/
    fmt::print(fg(fmt::color::white_smoke), "{}->Number of dim/s:\t{}\n", party_name, (double(len_hist) / t) * 1e6);
    fmt::print(fg(fmt::color::white_smoke), "{}->Computing Time\t{} ms\n", party_name, t / (1000.0));
    fmt::print(fg(fmt::color::white_smoke), "{}->Computing Bytes Sent\t{} bytes, {} MB\n", party_name, total_comm, total_comm / (1024.0 * 1024.0));
    logger->info("{}->Number of dim/s:\t{}", party_name, (double(len_hist) / t) * 1e6);
    logger->info("{}->Computing Time\t{} ms", party_name, t / (1000.0));
    logger->info("{}->Computing Bytes Sent\t{} bytes, {} MB", party_name, total_comm, total_comm / (1024.0 * 1024.0));
#endif

    /******************* Cleanup ****************/
    spdlog::drop_all();
    return sharing_hist_vec_star;
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
    // int bw_compare = -1;
    // double percentile = -1.0;
    uint64_t k_double_prime = -1;
    int l_double_prime = -1;

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
    // amap.arg("bc", bw_compare, "Bitwidth of compare input");
    amap.arg("min", generate_min, "生成的数据的最小值");
    amap.arg("max", generate_max, "生成的数据的最大值");
    amap.arg("nd", number_of_data, "生成的数据的个数");
    amap.arg("kr", k_double_prime, "k_double_prime");
    amap.arg("lr", l_double_prime, "l_double_prime");
    amap.parse(argc, argv);

    my::set_value_int(len_hist, 1000);
    my::set_value_int(left, 100);
    my::set_value_int(step, 1);

    my::set_value_int(number_of_data, 500);
    my::set_value_uint64_t(generate_min, left + len_hist * step * 0.02);
    my::set_value_uint64_t(generate_max, left + (len_hist - 1) * step - len_hist * step * 0.02);

    // my::set_value_double(percentile, 0.1);
    my::set_value_uint64_t(k_double_prime, 1);
    my::set_value_int(l_double_prime, 2);

    int bw_hist_part1_1 = ceil(log2(number_of_data));
    int bw_hist_part1_2 = ceil(log2((len_hist - 1) * step)+1);
    int bw_hist_part1 = bw_hist_part1_1 > bw_hist_part1_2 ? bw_hist_part1_1 : bw_hist_part1_2;
    int bw_hist_part2 = ceil(log2(k_double_prime)+1) > l_double_prime ? ceil(log2(k_double_prime)+1) : l_double_prime;

    int bw_hist_1 = bw_hist_part1 + bw_hist_part2 + 2 + 1;
    int bw_hist_2 = ceil(log2(left + (len_hist - 1) * step)) + 2;

    my::set_value_int(bw_hist, bw_hist_1 > bw_hist_2 ? bw_hist_1 : bw_hist_2);

    my::set_value_int(bw_numeric, ceil(log2(left + (len_hist - 1) * step)));
    // my::set_value_int(bw_compare, 14);

    // my::set_value_int(bw_hist, 16);
    // my::set_value_int(bw_numeric, 16);
    // my::set_value_int(bw_compare, 14);

    log_filename = fmt::format("logs/my_ring_outlier_{}.txt", party);

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
    auto qinfo_1_4 = quantile_number(1, 2, number_of_data);
    auto qinfo_3_4 = quantile_number(3, 2, number_of_data);
    fmt::print(
        fg(fmt::color::blue),
        "-------------- {}-> outlier ----------------------\nnumber_of_data:{}, bw_numeric:{}, len_hist:{}, bw_hist:{}, left:{}, generate_min:{}, generate_max:{}, k_double_prime= {}, l_double_prime= {}, bw_hist_part1_1= {}, bw_hist_part1_2= {}, bw_hist_part2= {}\n",
        party_name, number_of_data, bw_numeric, len_hist, bw_hist, left, generate_min, generate_max, k_double_prime, l_double_prime, bw_hist_part1_1, bw_hist_part1_2, bw_hist_part2);
    logger->info(
        "-------------- {}-> outlier ----------------------\nnumber_of_data:{}, bw_numeric:{}, len_hist:{}, bw_hist:{}, left:{}, generate_min:{}, generate_max:{}, k_double_prime= {}, l_double_prime= {}, bw_hist_part1_1= {}, bw_hist_part1_2= {}, bw_hist_part2= {}\n",
        party_name, number_of_data, bw_numeric, len_hist, bw_hist, left, generate_min, generate_max, k_double_prime, l_double_prime, bw_hist_part1_1, bw_hist_part1_2, bw_hist_part2);

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

    std::vector<uint64_t> hist_vec_actual;
    if (party == ALICE)
    {
        auto numeric_data_vec = my::generate_numeric_data_vector(number_of_data, bw_numeric, generate_min, generate_max);
        //--------------------------------------------------------

#ifdef _VERIFY_
        vector<uint64_t> data_vec(numeric_data_vec);
        std::sort(data_vec.begin(), data_vec.end());
        quantile_1_4_actual.left = data_vec[qinfo_1_4.index - 1];
        quantile_1_4_actual.right = data_vec[qinfo_1_4.index];
        quantile_3_4_actual.left = data_vec[qinfo_3_4.index - 1];
        quantile_3_4_actual.right = data_vec[qinfo_3_4.index];
        fmt::print(fg(fmt::color::blue), "Line:{}, quantile_1_4_actual.left={}, quantile_1_4_actual.left={}, quantile_3_4_actual.left={}, quantile_3_4_actual.right={}\n", __LINE__, quantile_1_4_actual.left, quantile_1_4_actual.right, quantile_3_4_actual.left, quantile_3_4_actual.right);
#endif
        //--------------------------------------------------------

        hist_vec_actual = my::numeric_data_to_hist(len_hist, left, step, numeric_data_vec, number_of_data);
        std::vector<uint64_t> data_accumulate_vec(len_hist);
        data_accumulate_vec[0] = hist_vec_actual[0];
        for (int i = 1; i < len_hist; i++)
        {
            data_accumulate_vec[i] = data_accumulate_vec[i - 1] + hist_vec_actual[i];
        }
        sharing_hist_vec = my::generate_and_send_sharing_vector_plaintext(iopackArr[0], hist_vec_actual, bw_hist);

#ifdef _CHECK_
        for (int i = 1; i < len_hist; i++)
        {
            if (i < 25 || i >= len_hist - 10)
            {
                std::cout << "Debug2: data_accumulate_vec" << i << ": " << data_accumulate_vec[i] << std::endl;
            }
        }
#endif
    }
    else if (party == BOB)
    {
        sharing_hist_vec = my::recv_sharing_vector_plaintext_vector(iopackArr[0], len_hist);
    }

    ShareA sharing_Q_1_4 = hist_quantile(1, 2, sharing_hist_vec, left, step, bw_hist);
    ShareA sharing_Q_3_4 = hist_quantile(3, 2, sharing_hist_vec, left, step, bw_hist);
    auto sharing_V = outlier_detect(sharing_hist_vec, left, step, sharing_Q_3_4, sharing_Q_1_4, bw_hist, 2, k_double_prime, l_double_prime);
    std::vector<uint64_t> sharing_hist_star = outlier_estimation(sharing_hist_vec, sharing_V, bw_hist);

/************** Verification ****************/
#ifdef _VERIFY_

    auto hist_vec = my::reconstruct_uint64_t(party, sci::ALICE, iopackArr, sharing_hist_vec, bw_hist);
    auto hist_star = my::reconstruct_uint64_t(party, sci::ALICE, iopackArr, sharing_hist_star, bw_hist);

    if (party == BOB)
    {
        iopackArr[0]->io->send_data(&sharing_Q_1_4, sizeof(ShareA));
        iopackArr[0]->io->send_data(&sharing_Q_3_4, sizeof(ShareA));
        // iopackArr[0]->io->send_data(sharing_V.data(), sizeof(ShareA) * sharing_V.size());
        iopackArr[0]->io->send_data(sharing_V.data(), sizeof(ShareB) * sharing_V.size());
        // fmt::print(fg(fmt::color::light_yellow), "{}\n", sharing_V);
    }
    else if (party == ALICE)
    {
        ShareA sharing_Q_1_40;
        ShareA sharing_Q_3_40;
        auto sharing_V0(sharing_V);
        iopackArr[0]->io->recv_data(&sharing_Q_1_40, sizeof(ShareA));
        iopackArr[0]->io->recv_data(&sharing_Q_3_40, sizeof(ShareA));
        iopackArr[0]->io->recv_data(sharing_V0.data(), sizeof(ShareB) * sharing_V.size());

        // TODO 修改5
        uint64_t Q_1_4 = sci::unsigned_val(sharing_Q_1_4 + sharing_Q_1_40, bw_hist);
        uint64_t Q_3_4 = sci::unsigned_val(sharing_Q_3_4 + sharing_Q_3_40, bw_hist);
        double Q_1_4_double = Q_1_4 / pow(2, 2);
        double Q_3_4_double = Q_3_4 / pow(2, 2);
        double Q_1_4_double_actual = (1 - qinfo_1_4.gamma) * quantile_1_4_actual.left + qinfo_1_4.gamma * quantile_1_4_actual.right;
        double Q_3_4_double_actual = (1 - qinfo_3_4.gamma) * quantile_3_4_actual.left + qinfo_3_4.gamma * quantile_3_4_actual.right;
        fmt::print(fg(fmt::color::blue), "Line:{}, Q_1_4= {}, Q_1_4_double_actual= {}, Q_3_4= {}, Q_1_4_double_actual= {}\n", __LINE__, Q_1_4, Q_1_4_double_actual, Q_3_4, Q_3_4_double_actual);
        uint64_t IQR = Q_3_4 - Q_1_4;
        uint64_t x_H_actual = (Q_3_4 << l_double_prime) + k_double_prime * IQR;
        uint64_t x_L_actual = (Q_1_4 << l_double_prime) - k_double_prime * IQR;
        x_L_actual = int64_t(x_L_actual) < 0 ? 0 : x_L_actual;
        fmt::print(fg(fmt::color::blue), "Line:{}, x_H_actual= {}, x_L_actual= {}\n", __LINE__, x_H_actual, x_L_actual);
        assert(Q_1_4_double == Q_1_4_double_actual);
        assert(Q_3_4_double == Q_3_4_double_actual);

        std::vector<ShareB> V(len_hist);
        std::transform(std::execution::par_unseq,
                       sharing_V.begin(), sharing_V.end(), sharing_V0.begin(),
                       V.begin(), [](auto e1, auto e2)
                       { return (e1 + e2) % 2; });

        std::vector<uint64_t> dataspace(len_hist);
        uint64_t start_shifted = left << (l_double_prime + 2);
        uint64_t step_shifted = step << (l_double_prime + 2);
        std::generate(dataspace.begin(), dataspace.end(), [start_shifted, step_shifted]()
                      { static uint64_t i = start_shifted-step_shifted; return i+=step_shifted; });

        std::vector<uint8_t> V_actual(len_hist);
        std::transform(std::execution::par_unseq,
                       dataspace.begin(), dataspace.end(),
                       V_actual.begin(), [x_L_actual, x_H_actual](auto e)
                       { return ((e < x_L_actual) || (e > x_H_actual)); });
        std::vector<uint64_t> hist_star_actual(sharing_hist_vec.size());
        std::transform(std::execution::par_unseq,
                       V_actual.begin(), V_actual.end(), hist_vec_actual.begin(),
                       hist_star_actual.begin(), [](uint8_t sel, auto e)
                       { return sel ? 0 : e; });

#ifdef _CHECK_
        // fmt::print("hist_vec: {}\n", hist_vec);
        // fmt::print("hist_star: {}\n", hist_star);
        // fmt::print("hist_star_actual: {}\n", hist_star_actual);
        // fmt::print(fg(fmt::color::crimson), "Line:{}, V= {}\n", __LINE__, V);
        // fmt::print(fg(fmt::color::crimson), "Line:{}, V_actual= {}\n", __LINE__, V_actual);
#endif

        for (int i = 0; i < len_hist; i++)
        {
            assert(V_actual[i] == V[i]);
            assert(hist_star_actual[i] == hist_star[i]);
        }
    }
#endif

    for (int i = 0; i < num_threads; i++)
    {
        delete iopackArr[i];
        delete otpackArr[i];
    }
    spdlog::drop_all();
}
