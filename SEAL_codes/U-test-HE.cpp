#include "examples.h"
#include <bits/stdc++.h>
#include <random>
#include "ArgMapping.h"
#include <fmt/core.h>
#include <fmt/color.h>
#include <sstream>
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/async.h"
// #include "vector"

spdlog::filename_t log_filename;

// the histgram is randomly generated
std::vector<uint64_t> generateRandomUint64Array(int number_of_data, int modulus)
{
    std::vector<uint64_t> result;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> dis(0, modulus - 1);
    for (int i = 0; i < number_of_data; i++)
    {
        result.push_back(dis(gen));
    }
    return result;
}

int main(int argc, char **argv)
{
    int len_hist = 4000;
    int bw_hist = -1;
    int number_of_data = 400;
    size_t poly_modulus_degree = 8192;

    ArgMapping amap;
    amap.arg("lh", len_hist, "Length of histgram vector");
    amap.arg("bh", bw_hist, "Bitwidth of histgram");
    amap.arg("nd", number_of_data, "生成的数据的个数");
    amap.arg("n", poly_modulus_degree, "poly_modulus_degree");
    amap.parse(argc, argv);

    if (bw_hist == -1)
    {
        bw_hist = ceil(2*log2(number_of_data)+log2(3*len_hist));
    }
    fmt::print(fg(fmt::color::green), "bw_hist= {}\n", bw_hist);
    uint64_t plaintext_modulus = 1ULL << bw_hist;

    log_filename = fmt::format("logs/my_HE_U.txt");

    seal::EncryptionParameters parms(seal::scheme_type::bfv);
    parms.set_poly_modulus_degree(poly_modulus_degree);

    parms.set_coeff_modulus(seal::CoeffModulus::BFVDefault(poly_modulus_degree));
    parms.set_plain_modulus(plaintext_modulus);
    seal::SEALContext context(parms);
    print_parameters(context);
    seal::KeyGenerator keygen(context);
    
    double total_time = 0;
    int rounds = 20;
    for (size_t i = 0; i < rounds; i++)
    {

        seal::SecretKey secret_key = keygen.secret_key();
        seal::PublicKey public_key;
        keygen.create_public_key(public_key);

        seal::Encryptor encryptor(context, public_key);
        seal::Decryptor decryptor(context, secret_key);
        seal::Evaluator evaluator(context);
        // seal::BatchEncoder batch_encoder(context);

        auto HistX = generateRandomUint64Array(len_hist, plaintext_modulus);
        auto HistY = generateRandomUint64Array(len_hist, plaintext_modulus);

        seal::Plaintext HistX_pt = seal::Plaintext();
        HistX_pt.resize(poly_modulus_degree);
        seal::Plaintext HistY_pt = seal::Plaintext();
        HistY_pt.resize(poly_modulus_degree);

        std::copy(HistX.begin(), HistX.end(), HistX_pt.data_span().begin());
        std::copy(HistY.begin(), HistY.end(), HistY_pt.data_span().begin());

        seal::Ciphertext HistX_ct = seal::Ciphertext();
        seal::Ciphertext HistY_ct = seal::Ciphertext();
        encryptor.encrypt(HistX_pt, HistX_ct);
        encryptor.encrypt(HistY_pt, HistY_ct);

        using namespace std::chrono;
        auto evaluate_start_time = system_clock::now();

        std::vector<uint64_t> X_mult(HistX.size(), 2);
        std::vector<uint64_t> Y_mult(HistY.size(), 2);
        X_mult[X_mult.size() - 1] = 1;
        Y_mult[0] = 1;

        seal::Plaintext X_mult_pt;
        seal::Plaintext Y_mult_pt;
        X_mult_pt.resize(poly_modulus_degree);
        Y_mult_pt.resize(poly_modulus_degree);

        std::copy(X_mult.begin(), X_mult.end(), X_mult_pt.data_span().begin());
        std::copy(Y_mult.begin(), Y_mult.end(), Y_mult_pt.data_span().begin());

        seal::Ciphertext WX_ct;
        seal::Ciphertext WY_ct;

        evaluator.multiply_plain(HistX_ct, X_mult_pt, WX_ct);
        evaluator.multiply_plain(HistY_ct, Y_mult_pt, WY_ct);

        seal::Ciphertext UX_ct;
        seal::Ciphertext UY_ct;
        evaluator.multiply(WX_ct, HistY_ct, UX_ct);
        evaluator.multiply(WY_ct, HistX_ct, UY_ct);

        auto evaluate_end_time = system_clock::now();
        auto evaluate_duration = duration_cast<microseconds>(evaluate_end_time - evaluate_start_time);
        total_time += evaluate_duration.count();
    }

    fmt::print(fg(fmt::color::green), "evaluate_time= {} us\n", total_time / rounds);
}