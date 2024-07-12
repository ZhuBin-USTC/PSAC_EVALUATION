#include "Math/math-functions.h"
#include "MyFunctions/my-functions.h"
#include <iostream>
#include <thread>
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

// bw_hist1 = ceil(log2(number_of_data)) + l_prime + 1
// bw_hist2 = ceil(log2(left + (len_hist - 1) * step)) + l_prime
// bw_hist = max(bw_hist1 , bw_hist2)

/*
Q(p, vector x) = (1-gamma) * x_index + gamma * x_{index+1}
index = floor((n-1) * p) + 1
gamma = (n-1) * p - floor((n-1) * p)
n is the length of vector x
x begins from 1
*/
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
    amap.parse(argc, argv);

    my::set_value_int(len_hist, 1000);
    my::set_value_int(left, 100);
    my::set_value_int(step, 1);

    my::set_value_int(number_of_data, 500);
    my::set_value_uint64_t(generate_min, left + len_hist * step * 0.02);
    my::set_value_uint64_t(generate_max, left + (len_hist - 1) * step - len_hist * step * 0.02);

    // my::set_value_double(percentile, 0.1);
    my::set_value_uint64_t(p_prime, 1);
    my::set_value_int(l_prime, 2);
    int bw_hist1 = ceil(log2(number_of_data)) + l_prime + 1;
    int bw_hist2 = ceil(log2(left + (len_hist - 1) * step)) + l_prime;
    my::set_value_int(bw_hist, bw_hist1 > bw_hist2 ? bw_hist1 : bw_hist2);
    my::set_value_int(bw_numeric, ceil(log2(left + (len_hist - 1) * step)));

    log_filename = fmt::format("logs/my_ring_quantile_hideK_{}.txt", party);

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
               "-------------- {}-> quantile_hideK ----------------------\nnumber_of_data:{}, bw_numeric:{}, len_hist:{}, bw_hist:{}, left:{}, generate_min:{}, generate_max:{}, percentile:{}, p_prime= {}, l_prime= {}, bw_hist1= {},bw_hist2= {}\n",
               party_name, number_of_data, bw_numeric, len_hist, bw_hist, left, generate_min, generate_max, qinfo.percentile, p_prime, l_prime, bw_hist1, bw_hist2);
    logger->info(
        "-------------- {}-> quantile_hideK ----------------------\nnumber_of_data:{}, bw_numeric:{}, len_hist:{}, bw_hist:{}, left:{}, generate_min:{}, generate_max:{}, percentile:{}, p_prime= {}, l_prime= {}, bw_hist1= {},bw_hist2= {}\n",
        party_name, number_of_data, bw_numeric, len_hist, bw_hist, left, generate_min, generate_max, qinfo.percentile, p_prime, l_prime, bw_hist1, bw_hist2);

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
        fmt::print(fg(fmt::color::red), "Line:{}\n", __LINE__);

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

    for (int i = 0; i < num_threads; i++)
    {
        delete iopackArr[i];
        delete otpackArr[i];
    }
    spdlog::drop_all();
}
