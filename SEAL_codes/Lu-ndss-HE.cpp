#include "examples.h"
#include <bits/stdc++.h>
#include "Lu-ndss-HE.h"
#include <random>
#include "ArgMapping.h"
#include <fmt/core.h>
#include <fmt/color.h>
#include <sstream>

// #define PRINT_DBG
// #define TIME_COUNT

#define varName(x) #x

using namespace std;
using namespace seal;

NDSS::NDSS(seal::scheme_type scheme, size_t poly_modulus_degree, int bit_size) : parms(scheme), poly_modulus_degree(poly_modulus_degree)
{
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(seal::CoeffModulus::BFVDefault(poly_modulus_degree));
    parms.set_plain_modulus(seal::PlainModulus::Batching(poly_modulus_degree, bit_size));
    context = std::make_unique<seal::SEALContext>(parms);
    fmt::print(fg(fmt::color::green), "Set encryption parameters and print\n");
    print_parameters(*context.get());
    keygen = std::make_unique<seal::KeyGenerator>(*context.get());
    secret_key = keygen->secret_key();
    keygen->create_public_key(public_key);
    keygen->create_relin_keys(relin_keys);
    keygen->create_galois_keys(galois_keys);
    encryptor = std::make_unique<seal::Encryptor>(*context.get(), public_key);
    decryptor = std::make_unique<seal::Decryptor>(*context.get(), secret_key);
    evaluator = std::make_unique<seal::Evaluator>(*context.get());
    batch_encoder = std::make_unique<seal::BatchEncoder>(*context.get());
}

void NDSS::Repeat(
    Ciphertext &u, Ciphertext &u_hat, int theta, int R, std::size_t print_size) const
{
    int rho = ceil(log2(R + 1)) - 1;

#ifdef PRINT_DBG
    cout << "rho" << rho << endl;
#endif

    int k = theta;
    encryptor->encrypt_zero(u_hat);

    size_t slot_count = batch_encoder->slot_count();
    size_t row_size = slot_count / 2;

    vector<uint64_t> pod_matrix(slot_count, 0ULL);
    for (int i = 0; i < theta; i++)
    {
        pod_matrix[i] = 1ULL;
        pod_matrix[i + row_size] = 1ULL;
    }
    Plaintext x_plain;
    batch_encoder->encode(pod_matrix, x_plain);

#ifdef PRINT_DBG
    print_line(__LINE__);
    print_plaintext(batch_encoder, x_plain);
#endif

    Ciphertext u_copy;
    evaluator->multiply_plain(u, x_plain, u_copy);

#ifdef PRINT_DBG
    print_line(__LINE__);
    print_ciphertext(batch_encoder, decryptor, u_copy, print_size);
#endif

    for (int i = 0; i <= rho; i++)
    {
        if ((R % 2) == 1)
        {
#ifdef PRINT_DBG
            cout << "1";
#endif
            evaluator->rotate_rows_inplace(u_hat, -k, galois_keys);
            evaluator->add_inplace(u_hat, u_copy);
        }

#ifdef PRINT_DBG
        else
        {
            cout << "0";
        }
#endif

        Ciphertext rotated;
        evaluator->rotate_rows(u_copy, -k, galois_keys, rotated);
        evaluator->add(u_copy, rotated, u_copy);
        k = k + k;
        R = R >> 1;
    }
}

void repeat_plain(vector<uint64_t> &u, int theta, int R, vector<uint64_t> &u_hat)
{
    vector<uint64_t> u_copy = u;
    int rho = ceil(log2(R + 1)) - 1;
    cout << "rho" << rho << endl;
    int k = theta;
    for (int i = 0; i < u_copy.size(); i++)
    {
        u_hat[i] = 0;
    }
    for (int i = theta; i < u_copy.size(); i++)
    {
        u_copy[i] = 0;
    }
    int R_temp = R;
    for (int i = 0; i <= rho; i++)
    {
        if ((R_temp % 2) == 1)
        {
            cout << "1";
            vector_right_rotate(u_hat, k);
            for (int j = 0; j < u_copy.size(); j++)
            {
                u_hat[j] = u_hat[j] + u_copy[j];
            }
        }
        else
        {
            cout << "0";
        }
        auto temp = u_copy;
        vector_right_rotate(temp, k);
        for (int j = 0; j < u_copy.size(); j++)
        {
            u_copy[j] = u_copy[j] + temp[j];
        }
        k = k + k;
        R_temp = R_temp >> 1;
    }
}

std::vector<std::vector<size_t>> generate_random_permutations(int theta, int D)
{
    std::vector<std::vector<size_t>> pi(theta);
    pi[0] = std::vector<size_t>(D);
    for (int j = 0; j < D; j++)
    {
        pi[0][j] = j;
    }
    for (int i = 1; i < theta; i++)
    {
        pi[i] = pi[0];
    }
    for (int i = 0; i < theta; i++)
    {
        random_shuffle_vector(pi[i]);
    }
    return pi;
}

std::vector<uint64_t> generate_w_alpha_j(int theta, int D)
{
    std::vector<uint64_t> w_alpha_j(theta * D);
    auto pi = generate_random_permutations(theta, D);
    for (size_t alpha = 0; alpha < D; alpha++)
    {
        for (size_t j = 0; j < theta; j++)
        {
            w_alpha_j[theta * alpha + j] = pi[j][alpha];
        }
    }
    return w_alpha_j;
}

// 将向量array分成两半，放在matrix的两行中
std::vector<uint64_t> NDSS::array_to_SIMD_matrix(vector<uint64_t> &array)
{
    size_t array_size = array.size();
    size_t half_size = (array_size + 1) / 2;

    size_t slot_count = batch_encoder->slot_count();
    size_t row_size = slot_count / 2;
    std::vector<uint64_t> pod_matrix(slot_count, 0ULL);
    std::copy(array.begin(), array.begin() + half_size, pod_matrix.begin());
    std::copy(array.begin() + half_size, array.end(), pod_matrix.begin() + row_size);
    return pod_matrix;
}

seal::Plaintext NDSS::array_to_plaintext(vector<uint64_t> &array)
{
    auto pod_matrix = array_to_SIMD_matrix(array);
    Plaintext pt;
    batch_encoder->encode(pod_matrix, pt);
    return pt;
}

// 将输入向量分成两行放进SIMD编码的两行矩阵中
void NDSS::batch_greater_than(seal::Ciphertext &a_ct, seal::Ciphertext &b_ct, seal::Ciphertext &gamma_ct, const BGTContext &bGTContext, std::size_t print_size)
{
#ifdef TIME_COUNT
    using namespace std::chrono;
    auto start = system_clock::now();
#endif

    srand(time(NULL));
    seal::Ciphertext a_repeated;
    seal::Ciphertext b_repeated;
    int half_theta = bGTContext.half_theta;
    int D = bGTContext.D;
    size_t total_size = half_theta * D * 2;
    size_t half_size = half_theta * D;
    Repeat(a_ct, a_repeated, half_theta, D, print_size);
    Repeat(b_ct, b_repeated, half_theta, D, print_size);

    auto w_alpha_j0 = generate_w_alpha_j(half_theta, D);
    auto w_alpha_j1 = generate_w_alpha_j(half_theta, D);

    size_t slot_count = batch_encoder->slot_count();
    size_t row_size = slot_count / 2;
    vector<uint64_t> w_pod_matrix(slot_count, 0ULL);

#ifdef PRINT_DBG
    print_line(__LINE__);
    print_vector(w_alpha_j0, print_size);
    print_vector(w_alpha_j1, print_size);
#endif
    std::copy(w_alpha_j0.begin(), w_alpha_j0.end(), w_pod_matrix.begin());
    std::copy(w_alpha_j1.begin(), w_alpha_j1.end(), w_pod_matrix.begin() + row_size);

    seal::Plaintext w_alpha_j_pt;
    batch_encoder->encode(w_pod_matrix, w_alpha_j_pt);

    seal::Ciphertext beta_ct;
    evaluator->sub(a_repeated, b_repeated, beta_ct);

#ifdef TIME_COUNT
    auto finish_rotate = system_clock::now();
    auto duration_rotate = duration_cast<microseconds>(finish_rotate - start);
    auto cost_rotate = double(duration_rotate.count()) * microseconds::period::num / microseconds::period::den;
    fmt::print(fg(fmt::color::green), "time rotate: {} s\n", cost_rotate);
    auto start2 = system_clock::now();
#endif

    evaluator->sub_plain_inplace(beta_ct, w_alpha_j_pt);
#ifdef PRINT_DBG
    print_line(__LINE__);
    print_plaintext(w_alpha_j_pt, print_size);
    print_ciphertext(beta_ct, print_size);
#endif
    uint64_t T = parms.plain_modulus().value();
    vector<uint64_t> r(total_size, 0ULL);

    std::default_random_engine e;
    std::uniform_int_distribution<uint64_t> u(1, T - 1);
    e.seed(time(0));

    for (int i = 0; i < total_size; i++)
    {
        r[i] = u(e);
    }
#ifdef PRINT_DBG
    print_line(__LINE__);
    print_vector(r, "r");
#endif
    seal::Plaintext r_pt = array_to_plaintext(r);
    // print_line(__LINE__);
    // print_plaintext(r_pt);
    evaluator->multiply_plain(beta_ct, r_pt, gamma_ct);
#ifdef PRINT_DBG
    print_line(__LINE__);
    print_ciphertext(gamma_ct, print_size);
#endif
}

void NDSS::batch_greater_than_evaluate(seal::Ciphertext &gamma_ct, std::vector<bool> &result, const BGTContext &bGTContext, std::size_t print_size)
{
    srand(time(NULL));
    // int half_theta = bGTContext.half_theta;
    int D = bGTContext.D;
    size_t total_size = bGTContext.half_theta * D * 2;
    size_t half_size = bGTContext.half_theta * D;

    Plaintext gamma_pt;
    decryptor->decrypt(gamma_ct, gamma_pt);

    // size_t total_size = theta * D;
    size_t slot_count = batch_encoder->slot_count();
    size_t row_size = slot_count / 2;
    vector<uint64_t> gamma_pod_matrix(slot_count, 0ULL);
    batch_encoder->decode(gamma_pt, gamma_pod_matrix);
    // print_line(__LINE__);
    // print_matrix(gamma_pod_matrix, row_size, print_size);
    for (int i = 0; i < bGTContext.theta; i++)
    {
        result[i] = false;
    }

    // for (int i = 0; i < total_size; i++)
    // {
    //     if (gamma_pod_matrix[i] == 0)
    //     {
    //         result[i % theta] = true;
    //     }
    // }

    for (int i = 0; i < half_size; i++)
    {
        if (gamma_pod_matrix[i] == 0)
        {
            result[i % bGTContext.half_theta] = true;
        }
    }
    for (int i = row_size; i < row_size + half_size; i++)
    {
        if (gamma_pod_matrix[i] == 0)
        {
            size_t result_index = (i - row_size + half_size) % bGTContext.half_theta + bGTContext.half_theta;
            if (result_index < bGTContext.theta)
            {
                result[result_index] = true;
            }
        }
    }
}

void NDSS::large_batch_greater_than_evaluate(seal::Ciphertext &a_ct, seal::Ciphertext &b_ct, std::vector<seal::Ciphertext> &gamma_cts, const BGTContext &bGTContext, std::size_t print_size)
{
#ifdef TIME_COUNT
    using namespace std::chrono;
    auto start = system_clock::now();
#endif
    seal::Ciphertext a_repeated;
    seal::Ciphertext b_repeated;

    Repeat(a_ct, a_repeated, bGTContext.half_theta, bGTContext.repeat_times_each_ct, print_size);
    Repeat(b_ct, b_repeated, bGTContext.half_theta, bGTContext.repeat_times_each_ct, print_size);

    seal::Ciphertext beta_ct;
    evaluator->sub(a_repeated, b_repeated, beta_ct);

    auto w_alpha_j0 = generate_w_alpha_j(bGTContext.half_theta, bGTContext.D);
    auto w_alpha_j1 = generate_w_alpha_j(bGTContext.half_theta, bGTContext.D);

    std::vector<std::vector<uint64_t>> w_pod_matrixs(bGTContext.number_of_cts, std::vector<uint64_t>(bGTContext.slot_count, 0ULL));

#ifdef PRINT_DBG
    print_line(__LINE__);
    print_vector(w_alpha_j0, print_size);
    print_vector(w_alpha_j1, print_size);
#endif
    auto w0_begin = w_alpha_j0.begin();
    auto w1_begin = w_alpha_j1.begin();
    for (size_t i = 0; i < bGTContext.number_of_cts - 1; i++)
    {
        std::copy(w0_begin, w0_begin + bGTContext.used_each_row, w_pod_matrixs[i].begin());
        std::copy(w1_begin, w1_begin + bGTContext.used_each_row, w_pod_matrixs[i].begin() + bGTContext.row_size);
        w0_begin = w0_begin + bGTContext.used_each_row;
        w1_begin = w1_begin + bGTContext.used_each_row;
    }
    std::copy(w0_begin, w_alpha_j0.end(), w_pod_matrixs[bGTContext.number_of_cts - 1].begin());
    std::copy(w1_begin, w_alpha_j1.end(), w_pod_matrixs[bGTContext.number_of_cts - 1].begin() + bGTContext.row_size);

    std::vector<seal::Plaintext> w_alpha_j_pts(bGTContext.number_of_cts);

    for (int i = 0; i < bGTContext.number_of_cts; i++)
    {
        batch_encoder->encode(w_pod_matrixs[i], w_alpha_j_pts[i]);
    }

#ifdef PRINT_DBG
    for (int i = 0; i < bGTContext.number_of_cts; i++)
    {
        std::cout << "w_pod_pts[" << i << "]:";
        print_plaintext(w_alpha_j_pts[i], 100);
    }
#endif
    std::vector<seal::Ciphertext> beta_cts(bGTContext.number_of_cts, beta_ct);

#ifdef TIME_COUNT
    auto finish_rotate = system_clock::now();
    auto duration_rotate = duration_cast<microseconds>(finish_rotate - start);
    auto cost_rotate = double(duration_rotate.count()) * microseconds::period::num / microseconds::period::den;
    std::cout << "time rotate: " << cost_rotate << " s" << std::endl;

    auto start2 = system_clock::now();
#endif

    for (int i = 0; i < bGTContext.number_of_cts; i++)
    {
        evaluator->sub_plain_inplace(beta_cts[i], w_alpha_j_pts[i]);
    }

    std::vector<uint64_t> zeros(bGTContext.used_each_ct, 0ULL);
    std::vector<std::vector<uint64_t>> random_numbers(bGTContext.number_of_cts, zeros);
    std::vector<seal::Plaintext> random_pts(bGTContext.number_of_cts);

    std::default_random_engine e;
    uint64_t T = parms.plain_modulus().value();
    std::uniform_int_distribution<uint64_t> u(1, T - 1);
    e.seed(time(0));

    for (size_t i = 0; i < bGTContext.number_of_cts; i++)
    {
        for (size_t j = 0; j < bGTContext.used_each_ct; j++)
        {
            random_numbers[i][j] = u(e);
        }
#ifdef PRINT_DBG
        print_line(__LINE__);
        std::cout << "random_numbers[" << i << "]: ";
        print_vector(random_numbers[i], " ");
#endif
        random_pts[i] = array_to_plaintext(random_numbers[i]);
        evaluator->multiply_plain(beta_cts[i], random_pts[i], gamma_cts[i]);
        // evaluator->relinearize_inplace(gamma_cts[i], relin_keys);
    }
}

std::vector<bool> NDSS::large_batch_greater_than_decrypt(std::vector<seal::Ciphertext> &gamma_cts, const BGTContext &bGTContext, std::size_t print_size)
{
    std::vector<bool> result(bGTContext.theta, false);
    std::vector<seal::Plaintext> gamma_pts(bGTContext.number_of_cts);
    std::vector<vector<uint64_t>> gamma_pod_matrixs(bGTContext.number_of_cts, vector<uint64_t>(bGTContext.slot_count, 0ULL));
    for (int i = 0; i < bGTContext.number_of_cts; i++)
    {
        decryptor->decrypt(gamma_cts[i], gamma_pts[i]);
        batch_encoder->decode(gamma_pts[i], gamma_pod_matrixs[i]);
    }
    for (size_t j = 0; j < bGTContext.number_of_cts - 1; j++)
    {
        auto gamma_pod_matrix = gamma_pod_matrixs[j];
        for (int i = 0; i < bGTContext.used_each_row; i++)
        {
            if (gamma_pod_matrix[i] == 0)
            {
                result[i % bGTContext.half_theta] = true;
            }
        }
        for (int i = bGTContext.row_size; i < bGTContext.row_size + bGTContext.used_each_row; i++)
        {
            if (gamma_pod_matrix[i] == 0)
            {
                size_t result_index = (i - bGTContext.row_size) % bGTContext.half_theta + bGTContext.half_theta;
                if (result_index < bGTContext.theta)
                {
                    result[result_index] = true;
                }
            }
        }
    }
    auto gamma_pod_matrix = gamma_pod_matrixs[bGTContext.number_of_cts - 1];
    for (int i = 0; i < bGTContext.used_last_row; i++)
    {
        if (gamma_pod_matrix[i] == 0)
        {
            result[i % bGTContext.half_theta] = true;
        }
    }
    for (int i = bGTContext.row_size; i < bGTContext.row_size + bGTContext.used_last_row; i++)
    {
        if (gamma_pod_matrix[i] == 0)
        {
            size_t result_index = (i - bGTContext.row_size) % bGTContext.half_theta + bGTContext.half_theta;
            if (result_index < bGTContext.theta)
            {
                result[result_index] = true;
            }
        }
    }
    return result;
}

void test_1()
{
    int poly_modulus_degree = 8192;
    int theta_Oj = 100;
    int N = 500;
    //   ndss_time(8192,100,500);
    using namespace seal;
    using namespace std::chrono;

    auto start0 = system_clock::now();

    NDSS ndss(scheme_type::bgv, poly_modulus_degree, 20);
    size_t slot_count = ndss.batch_encoder->slot_count();
    // int theta = theta_Oj;
    // int D = N;
    BGTContext bGTContext(slot_count, theta_Oj, N);

    // 模拟生成用户密文
    std::vector<Ciphertext> cts(N);
    std::default_random_engine e;
    std::uniform_int_distribution<size_t> u(0, theta_Oj - 1);
    e.seed(time(0));
    int i = 0;
    for (auto &ct : cts)
    {
        size_t indicator = u(e);
        std::vector<uint64_t> indicator_array(bGTContext.theta, 0ULL);
        std::fill(indicator_array.begin() + indicator, indicator_array.end(), 1);

        Plaintext pt = ndss.array_to_plaintext(indicator_array);
        ndss.encryptor->encrypt(pt, ct);
    }
    //---------------------------------------------------------

    auto finish0 = system_clock::now();
    auto duration0 = duration_cast<microseconds>(finish0 - start0);
    auto cost0 = double(duration0.count()) * microseconds::period::num / microseconds::period::den;
    std::cout << duration0.count() << " " << microseconds::period::num << " " << microseconds::period::den << "time: " << cost0 << " s" << std::endl;

    auto start = system_clock::now();

    Ciphertext ct_agg;
    ndss.encryptor->encrypt_zero(ct_agg);

    for (auto const &ct : cts)
    {
        ndss.evaluator->add_inplace(ct_agg, ct);
    }
#ifdef PRINT_DBG
    print_line(__LINE__);
    auto pt_agg = ndss.print_ciphertext(ct_agg, 100);
#endif
    int k = 30;
    std::vector<uint64_t> k_array(bGTContext.theta, ceil(double(k * N) / 100));

    Plaintext k_pt = ndss.array_to_plaintext(k_array);
    Ciphertext k_ct;

    ndss.encryptor->encrypt(k_pt, k_ct);

    std::vector<seal::Ciphertext> gamma_cts(bGTContext.number_of_cts);

    ndss.large_batch_greater_than_evaluate(ct_agg, k_ct, gamma_cts, bGTContext);
    auto result = ndss.large_batch_greater_than_decrypt(gamma_cts, bGTContext, 100);
    print_vector(result, "result:");
    for (int i = 1; i < bGTContext.theta; i++)
    {
        if (result[i] != result[i - 1])
            std::cout << i << std::endl;
    }
}

void ndss_time(int poly_modulus_degree, int theta_Oj, int N, int k)
{

    using namespace seal;
    using namespace std::chrono;
    int K = ceil(double(k * N) / 100);

    NDSS ndss(scheme_type::bgv, poly_modulus_degree, 20);
    size_t slot_count = ndss.batch_encoder->slot_count();
    BGTContext bGTContext(slot_count, theta_Oj, N);

    // 模拟生成用户密文
    std::vector<Ciphertext> cts(N);
    std::default_random_engine e;
    std::uniform_int_distribution<size_t> u(0, theta_Oj - 1);
    e.seed(time(0));

    /*
    int i = 0;
    for (auto &ct : cts)
    {
      size_t indicator = u(e);
      std::vector<uint64_t> indicator_array(bGTContext.theta, 0ULL);
      std::fill(indicator_array.begin() + indicator, indicator_array.end(), 1);

      Plaintext pt = ndss.array_to_plaintext(indicator_array);
      ndss.encryptor->encrypt(pt, ct);
    }

    auto agg_start_time = system_clock::now();

    Ciphertext ct_agg;
    ndss.encryptor->encrypt_zero(ct_agg);

    for (auto const &ct : cts)
    {
      ndss.evaluator->add_inplace(ct_agg, ct);
    }

    auto agg_end_time = system_clock::now();
    auto agg_duration = duration_cast<microseconds>(agg_end_time - agg_start_time);
    auto agg_cost = double(agg_duration.count()) * microseconds::period::num / microseconds::period::den;
    std::cout << "agg_time: " << agg_cost << " s" << std::endl;

    */

    //---------------------------------------------------------

    std::vector<uint64_t> pt_agg_array(bGTContext.theta, 0ULL);
    int step = N / bGTContext.theta;
    for (int i = 0; i < bGTContext.theta; i++)
    {
        pt_agg_array[i] = i * step;
    }
    pt_agg_array[bGTContext.theta - 1] = N;
    Plaintext pt_agg = ndss.array_to_plaintext(pt_agg_array);
    Ciphertext ct_agg;
    ndss.encryptor->encrypt(pt_agg, ct_agg);

    auto evaluate_start_time = system_clock::now();

    std::vector<uint64_t> k_array(bGTContext.theta, K);
    Plaintext k_pt = ndss.array_to_plaintext(k_array);
    Ciphertext k_ct;
    ndss.encryptor->encrypt(k_pt, k_ct);
    std::vector<seal::Ciphertext> gamma_cts(bGTContext.number_of_cts);

    ndss.large_batch_greater_than_evaluate(ct_agg, k_ct, gamma_cts, bGTContext);
    // evaluate time

    vector<seal_byte> gamma_cts_bytes(bGTContext.number_of_cts * gamma_cts[0].save_size(compr_mode_type::none));
    for (int i = 0; i < bGTContext.number_of_cts; i++)
    {
        gamma_cts[i].save(gamma_cts_bytes.data() + i * gamma_cts[0].save_size(compr_mode_type::none), gamma_cts[0].save_size(compr_mode_type::none));
    }
    auto evaluate_end_time = system_clock::now();
    auto evaluate_duration = duration_cast<microseconds>(evaluate_end_time - evaluate_start_time);

    fmt::print(fg(fmt::color::green), "msg send: {} MB\n",
               double(bGTContext.number_of_cts * gamma_cts[0].save_size(compr_mode_type::none)) / (1024 * 1024));

    auto decrypt_start_time = system_clock::now();
    std::vector<seal::Ciphertext> gamma_cts_recv(bGTContext.number_of_cts);
    for (int i = 0; i < bGTContext.number_of_cts; i++)
    {
        // gamma_cts[i].save(gamma_cts_bytes.data()+i*gamma_cts[0].save_size(), gamma_cts[0].save_size());
        gamma_cts_recv[i].load(*ndss.context, gamma_cts_bytes.data() + i * gamma_cts[0].save_size(compr_mode_type::none), gamma_cts[0].save_size(compr_mode_type::none));
    }
    auto result = ndss.large_batch_greater_than_decrypt(gamma_cts_recv, bGTContext, 100);

    auto decrypt_end_time = system_clock::now();
    auto decrypt_duration = duration_cast<microseconds>(decrypt_end_time - decrypt_start_time);


    fmt::print(fg(fmt::color::green), "num cts= {}, evaluate_time= {} us, decrypt_time= {} us\n",
               gamma_cts.size(), evaluate_duration.count(), decrypt_duration.count());


    fmt::print(fg(fmt::color::red),
               "gamma_cts[0].size: {}, gamma_cts[0].size_capaci: {}\n",
               gamma_cts[0].size(), gamma_cts[0].size_capacity());
    fmt::print(fg(fmt::color::red),
               "size(seal_byte) = {}, sizeof(std::uint64_t) = {}\n",
               sizeof(seal_byte), sizeof(std::uint64_t));
    fmt::print(fg(fmt::color::red), "gamma_cts[0].save_size(compr_mode_type::none): {} bytes\n", gamma_cts[0].save_size(compr_mode_type::none));
    fmt::print(fg(fmt::color::red), "gamma_cts[0].save_size(compr_mode_type::zstd)): {} bytes\n", gamma_cts[0].save_size(compr_mode_type::zstd));
    fmt::print(fg(fmt::color::red), "gamma_cts[0].save_size(compr_mode_type::zlib): {} bytes\n", gamma_cts[0].save_size(compr_mode_type::zlib));
    // auto evaluate_cost = double(evaluate_duration.count()) * microseconds::period::num / microseconds::period::den;
    // auto decrypt_cost = double(decrypt_duration.count()) * microseconds::period::num / microseconds::period::den;
}

int main(int argc, char **argv)
{
    int poly_modulus_degree = 8192;
    int theta_Oj = 4000;
    int N = 400;
    int k = 30;
    ArgMapping amap;
    amap.arg("p", poly_modulus_degree, "poly_modulus_degree");
    amap.arg("O", theta_Oj, "theta_Oj");
    amap.arg("N", N, "Number of workers");
    amap.arg("k", k, "k-percentile");
    amap.parse(argc, argv);
    fmt::print(fg(fmt::color::green), "poly_modulus_degree: {}, theta_Oj: {}, N: {}, k: {}\n", poly_modulus_degree, theta_Oj, N, k);
    ndss_time(poly_modulus_degree, theta_Oj, N, k);
}