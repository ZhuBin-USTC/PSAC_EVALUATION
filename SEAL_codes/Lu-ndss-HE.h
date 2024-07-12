#pragma once
#include "examples.h"
#include "chrono"

struct BGTContext{
    size_t slot_count;
    int theta;
    int D;
    size_t row_size;
    int repeat_times_each_ct;
    int number_of_cts;
    size_t half_theta;
    size_t used_each_row;
    int repeat_times_last_ct;
    size_t used_last_row;
    size_t used_each_ct;
    BGTContext(size_t slot_count, int theta, int D):slot_count(slot_count),theta(theta),D(D)
    {
        row_size = slot_count / 2;
        repeat_times_each_ct = slot_count / theta;
        number_of_cts = ceil(double(D) / repeat_times_each_ct);
        half_theta = (theta + 1) / 2;
        used_each_row = half_theta * repeat_times_each_ct;
        used_each_ct = used_each_row * 2;
        repeat_times_last_ct = D - (number_of_cts - 1) * repeat_times_each_ct;
        used_last_row = half_theta * repeat_times_last_ct;
    }
};

class NDSS
{
public:
    size_t poly_modulus_degree;
    seal::EncryptionParameters parms;
    std::unique_ptr<seal::SEALContext> context;
    std::unique_ptr<seal::KeyGenerator> keygen;
    seal::SecretKey secret_key;
    seal::PublicKey public_key;
    seal::RelinKeys relin_keys;
    seal::GaloisKeys galois_keys;
    std::unique_ptr<seal::Encryptor> encryptor;
    std::unique_ptr<seal::Evaluator> evaluator;
    std::unique_ptr<seal::Decryptor> decryptor;
    std::unique_ptr<seal::BatchEncoder> batch_encoder;

    NDSS(seal::scheme_type scheme, size_t poly_modulus_degree, int bit_size = 20);
    void Repeat(seal::Ciphertext &u, seal::Ciphertext &u_hat, int theta, int R, std::size_t print_size = 25) const;

    std::vector<uint64_t> array_to_SIMD_matrix(std::vector<uint64_t> &array);
    seal::Plaintext array_to_plaintext(std::vector<uint64_t> &array);

    void batch_greater_than(seal::Ciphertext &a_ct, seal::Ciphertext &b_ct, seal::Ciphertext &gamma_ct, const BGTContext &bGTContext, std::size_t print_size = 25);
    void batch_greater_than_evaluate(seal::Ciphertext &gamma_ct, std::vector<bool> &result,const BGTContext &bGTContext, std::size_t print_size = 25);

    void large_batch_greater_than_evaluate(seal::Ciphertext &a_ct, seal::Ciphertext &b_ct, std::vector<seal::Ciphertext> &gamma_cts, const BGTContext &bGTContext, std::size_t print_size = 25);
    std::vector<bool> large_batch_greater_than_decrypt(std::vector<seal::Ciphertext> &gamma_cts, const BGTContext &bGTContext, std::size_t print_size = 25);

    void print_plaintext(seal::Plaintext &pt, std::size_t print_size = 25)
    {
        size_t slot_count = batch_encoder->slot_count();
        size_t row_size = slot_count / 2;
        std::vector<uint64_t> pod_result;
        batch_encoder->decode(pt, pod_result);
        print_matrix(pod_result, row_size, print_size);
    }

    std::vector<uint64_t> print_ciphertext(seal::Ciphertext &ct, std::size_t print_size = 25)
    {
        size_t slot_count = batch_encoder->slot_count();
        size_t row_size = slot_count / 2;
        seal::Plaintext decrypted_result;
        decryptor->decrypt(ct, decrypted_result);
        std::vector<uint64_t> pod_result;
        batch_encoder->decode(decrypted_result, pod_result);
        print_matrix(pod_result, row_size, print_size);
        return pod_result;
    }
};

void repeat_plain(std::vector<uint64_t> &u, int theta, int R, std::vector<uint64_t> &u_hat);
// void Repeat(const NDSS &ndss, seal::Ciphertext &u, seal::Ciphertext &u_hat, int theta, int R,            std::size_t print_size = 25);

template <typename T>
inline void print_vector(std::vector<T> &array, std::string name)
{
    std::cout << name << " [";
    for (auto i : array)
    {
        std::cout << i << ", ";
    }
    std::cout << "] ---> size=" << array.size() << std::endl;
}

template <typename T>
inline void random_shuffle_vector(std::vector<T> &input)
{
    random_shuffle(input.begin(), input.end());
}

inline void print_plain_array(seal::Plaintext plain, std::string name)
{
    std::cout << name << " [";
    auto &dynarray = plain.dyn_array();
    for (auto i : dynarray)
    {
        std::cout << i << ", ";
    }
    std::cout << "] ---> size=" << dynarray.size() << std::endl;
}

inline void vector_right_rotate(std::vector<uint64_t> &input, int rotate_size)
{
    std::rotate(input.begin(), input.end() - rotate_size, input.end());
}

inline void vector_left_rotate(std::vector<uint64_t> &input, int rotate_size)
{
    std::rotate(input.begin(), input.begin() + rotate_size, input.end());
}

std::vector<std::vector<size_t>> generate_random_permutations(int theta, int D);
std::vector<uint64_t> generate_w_alpha_j(int theta, int D);

void ndss_time(int poly_modulus_degree, int Oj, int N);