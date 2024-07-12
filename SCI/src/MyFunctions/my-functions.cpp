#include "MyFunctions/my-functions.h"
#include "my-functions.h"

// #define _VERIFY_

void my::set_value_int(int &var, int value)
{
    if (var == -1)
    {
        var = value;
    }
}

void my::set_value_uint64_t(uint64_t &var, uint64_t value)
{
    if (var == -1)
    {
        var = value;
    }
}

void my::set_value_double(double &var, double value)
{
    if (var == -1.0)
    {
        var = value;
    }
}

void my::setup_io_and_ots(int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[], int num_threads, int port, std::string address)
{
    for (int i = 0; i < num_threads; i++)
    {
        iopackArr[i] = new sci::IOPack(party, port + i, address);
        if (i & 1)
        {
            otpackArr[i] = new sci::OTPack(iopackArr[i], 3 - party);
        }
        else
        {
            otpackArr[i] = new sci::OTPack(iopackArr[i], party);
        }
    }
};

std::vector<ShareA> my::append_vectors(const std::vector<ShareA> input1, const std::vector<ShareA> input2)
{
    std::vector<ShareA> output(input1.size() + input2.size());
    std::copy(input1.begin(), input1.end(), output.begin());
    std::copy(input2.begin(), input2.end(), output.begin() + input1.size());
    return output;
}

ShareA my::shareA_const(int party, uint64_t x)
{
    if (party == sci::ALICE)
    {
        return x;
    }
    else
    {
        return 0;
    }
}

ShareB my::shareB_const(int party, uint8_t x)
{
    if (party == sci::ALICE)
    {
        return x;
    }
    else
    {
        return 0;
    }
}

uint64_t my::mask_gen_uint64_t(int bw)
{
    return (bw == 64 ? -1 : ((1ULL << bw) - 1));
}

std::vector<uint64_t> my::reconstruct_uint64_t(int party, int reconstruct_party, sci::IOPack *iopackArr[], std::vector<uint64_t> &sharing, int bw)
{
    auto mask = my::mask_gen_uint64_t(bw);
    std::vector<uint64_t> sharing_0(sharing);
    if (party == (3 - reconstruct_party))
    {
        iopackArr[0]->io->send_data(sharing.data(), sizeof(uint64_t) * sharing.size());
    }
    else if (party == reconstruct_party)
    {
        iopackArr[0]->io->recv_data(sharing_0.data(), sizeof(uint64_t) * sharing_0.size());
        std::transform(std::execution::par_unseq, sharing.begin(), sharing.end(), sharing_0.begin(), sharing_0.begin(), [mask](auto e1, auto e2)
                       { return (e1 + e2) & mask; });
    }
    return sharing_0;
}

std::vector<uint8_t> my::reconstruct_uint8_t(int party, int reconstruct_party, sci::IOPack *iopackArr[], std::vector<uint8_t> &sharing)
{
    std::vector<uint8_t> sharing_0(sharing);
    if (party == (3 - reconstruct_party))
    {
        iopackArr[0]->io->send_data(sharing.data(), sizeof(uint8_t) * sharing.size());
    }
    else if (party == reconstruct_party)
    {
        iopackArr[0]->io->recv_data(sharing_0.data(), sizeof(uint8_t) * sharing_0.size());
        std::transform(std::execution::par_unseq, sharing.begin(), sharing.end(), sharing_0.begin(), sharing_0.begin(), std::bit_xor());
    }
    return sharing_0;
}

// 专门用于本地测试的，不加HE协议
std::unique_ptr<uint64_t[]> my::generate_and_send_sharing_vector_plaintext(sci::IOPack *iopack, std::unique_ptr<uint64_t[]> &data_vector_ptr, int dim, int bitlength)
{
    sci::PRG128 prg;
    auto sharing_vector_ptr = std::make_unique<uint64_t[]>(dim);
    auto sharing_vector_send_ptr = std::make_unique<uint64_t[]>(dim);
    uint64_t mask_data = (bitlength == 64 ? -1 : ((1ULL << bitlength) - 1));
    // 给BOB生成sharing
    prg.random_data(sharing_vector_send_ptr.get(), dim * sizeof(uint64_t));
    for (int i = 0; i < dim; i++)
    {
        sharing_vector_send_ptr[i] = sharing_vector_send_ptr[i] & mask_data;
        // 给自己生成sharing
        sharing_vector_ptr[i] = (data_vector_ptr[i] - sharing_vector_send_ptr[i]) & mask_data;

#ifdef _VERIFY_
        if (i < 25)
        {
            std::cout << "> i=" << i << " " << sharing_vector_ptr[i] << " + " << sharing_vector_send_ptr[i] << " = " << data_vector_ptr[i] << std::endl;
        }
#endif
    }
    iopack->io->send_data(sharing_vector_send_ptr.get(), dim * sizeof(uint64_t));
    return sharing_vector_ptr;
}

std::vector<ShareA> my::generate_and_send_sharing_vector_plaintext(sci::IOPack *iopack, std::vector<uint64_t> &data_vector, int bitlength)
{
    sci::PRG128 prg;
    size_t dim = data_vector.size();
    auto sharing_vector = std::vector<ShareA>(dim);
    auto sharing_send_vector = std::vector<ShareA>(dim);

    uint64_t mask_data = mask_gen_uint64_t(bitlength);
    // 给BOB生成sharing
    prg.random_data(sharing_send_vector.data(), dim * sizeof(ShareA));

    // 给自己生成sharing
    std::transform(sharing_send_vector.begin(), sharing_send_vector.end(), sharing_send_vector.begin(),
                   [mask_data](ShareA sharing_send_element) -> ShareA
                   { return (sharing_send_element & mask_data); });

    std::transform(data_vector.begin(), data_vector.end(), sharing_send_vector.begin(), sharing_vector.begin(),
                   [mask_data](ShareA data_element, ShareA sharing_send_element) -> ShareA
                   { return ((data_element - sharing_send_element) & mask_data); });
    for (int i = 0; i < dim; i++)
    {

#ifdef _VERIFY_
        if (i < 25)
        {
            std::cout << "> i=" << i << " " << sharing_vector_ptr[i] << " + " << sharing_vector_send_ptr[i] << " = " << data_vector_ptr[i] << std::endl;
        }
#endif
    }
    iopack->io->send_data(sharing_send_vector.data(), dim * sizeof(uint64_t));
    return sharing_vector;
}

std::unique_ptr<uint64_t[]> my::recv_sharing_vector_plaintext(sci::IOPack *iopack, int dim)
{
    auto sharing_vector_ptr = std::make_unique<uint64_t[]>(dim);
    iopack->io->recv_data(sharing_vector_ptr.get(), dim * sizeof(uint64_t));
    return sharing_vector_ptr;
}

std::vector<uint64_t> my::recv_sharing_vector_plaintext_vector(sci::IOPack *iopack, int dim)
{
    std::vector<uint64_t> sharing_vector(dim);
    iopack->io->recv_data(sharing_vector.data(), dim * sizeof(uint64_t));
    return sharing_vector;
}

std::unique_ptr<uint64_t[]> my::generated_worker_data_arithmetic(int num_of_data, int bitlength_arithmetic, uint64_t min, uint64_t max)
{
    sci::PRG128 prg;
    std::unique_ptr<uint64_t[]> data_arithmetic_ptr = std::make_unique<uint64_t[]>(num_of_data);
    prg.random_data(data_arithmetic_ptr.get(), num_of_data * sizeof(uint64_t));
    uint64_t mask_arithmetic = (bitlength_arithmetic == 64 ? -1 : ((1ULL << bitlength_arithmetic) - 1));
    for (int i = 0; i < num_of_data; i++)
    {
        data_arithmetic_ptr[i] = (data_arithmetic_ptr[i] % (max - min + 1)) + min;
        data_arithmetic_ptr[i] = data_arithmetic_ptr[i] & mask_arithmetic;
    }
    return data_arithmetic_ptr;
}

std::vector<uint64_t> my::generate_numeric_data_vector(int num_of_data, int bitlength_arithmetic, uint64_t min, uint64_t max)
{
    sci::PRG128 prg;
    std::vector<uint64_t> data_arithmetic(num_of_data);
    prg.random_data(data_arithmetic.data(), num_of_data * sizeof(uint64_t));

    uint64_t mask_arithmetic = mask_gen_uint64_t(bitlength_arithmetic);

    std::transform(data_arithmetic.begin(), data_arithmetic.end(), data_arithmetic.begin(),
                   [max, min, mask_arithmetic](uint64_t e) -> uint64_t
                   { return ((e % (max - min + 1)) + min) & mask_arithmetic; });
    // for (int i = 0; i < num_of_data; i++)
    // {
    //     data_arithmetic[i] = (data_arithmetic[i] % (max - min + 1)) + min;
    //     data_arithmetic[i] = data_arithmetic[i] & mask_arithmetic;
    // }
    return data_arithmetic;
}

std::unique_ptr<uint64_t[]> my::data_arithmetic_to_frequency(
    int dim_frequency,
    uint64_t left, uint64_t step,
    std::unique_ptr<uint64_t[]> &data_arithmetic_ptr,
    int num_of_data)
{
    std::unique_ptr<uint64_t[]> data_frequency_ptr = std::make_unique<uint64_t[]>(dim_frequency);
    memset(data_frequency_ptr.get(), 0, sizeof(uint64_t) * dim_frequency);
    for (int i = 0; i < num_of_data; i++)
    {
        int index = (data_arithmetic_ptr[i] - left) / step;
        data_frequency_ptr[index]++;
    }
    return data_frequency_ptr;
}

std::vector<uint64_t> my::numeric_data_to_hist(
    int dim_frequency,
    uint64_t left, uint64_t step,
    std::vector<uint64_t> &data_arithmetic,
    int num_of_data)
{
    std::vector<uint64_t> data_frequency(dim_frequency, 0);
    // memset(data_frequency.data(), 0, sizeof(uint64_t) * dim_frequency);
    for (int i = 0; i < num_of_data; i++)
    {
        int index = (data_arithmetic[i] - left) / step;
        data_frequency[index]++;
    }
    return data_frequency;
}

void my::equal_thread(int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[],
                      int tid, uint8_t *res_eq, uint64_t *data, uint64_t value, int num_ops, int bitlength)
{
    std::unique_ptr<MyFunctions> myfunction_ptr;
    if (tid & 1)
    {
        myfunction_ptr = std::make_unique<MyFunctions>(3 - party, iopackArr[tid], otpackArr[tid]);
    }
    else
    {
        myfunction_ptr = std::make_unique<MyFunctions>(party, iopackArr[tid], otpackArr[tid]);
    }
    myfunction_ptr->check_equality(res_eq, data, value, num_ops, bitlength);
}

void my::equal_arithmetic_thread(int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[],
                                 int tid, uint8_t *res_eq, uint64_t *res_eq_arithmetic, uint64_t *data, uint64_t value, int num_ops, int bitlength_data, int bitlength_res)
{
    std::unique_ptr<MyFunctions> myfunction_ptr;
    if (tid & 1)
    {
        myfunction_ptr = std::make_unique<MyFunctions>(3 - party, iopackArr[tid], otpackArr[tid]);
    }
    else
    {
        myfunction_ptr = std::make_unique<MyFunctions>(party, iopackArr[tid], otpackArr[tid]);
    }
    myfunction_ptr->check_equality_arithmetic(res_eq, res_eq_arithmetic, data, value, num_ops, bitlength_data, bitlength_res);
}

void my::compare_const_thread(int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[],
                              int tid, uint8_t *res_cmp, uint64_t *data, uint64_t value, int num_ops, int bitlength, bool greater_than, bool equality)
{
    std::unique_ptr<MyFunctions> myfunction_ptr;
    if (tid & 1)
    {
        myfunction_ptr = std::make_unique<MyFunctions>(3 - party, iopackArr[tid], otpackArr[tid]);
    }
    else
    {
        myfunction_ptr = std::make_unique<MyFunctions>(party, iopackArr[tid], otpackArr[tid]);
    }
    myfunction_ptr->compare_const(res_cmp, data, value, num_ops, bitlength, greater_than, equality);
}

void my::compare_const_arithmetic_thread(int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[],
                                         int tid, uint8_t *res_cmp, uint64_t *res_cmp_arithmetic, uint64_t *data, uint64_t value, int num_ops, int bitlength_data, int bitlength_res, bool greater_than, bool equality)
{
    std::unique_ptr<MyFunctions> myfunction_ptr;
    if (tid & 1)
    {
        myfunction_ptr = std::make_unique<MyFunctions>(3 - party, iopackArr[tid], otpackArr[tid]);
    }
    else
    {
        myfunction_ptr = std::make_unique<MyFunctions>(party, iopackArr[tid], otpackArr[tid]);
    }
    myfunction_ptr->compare_const_arithmetic(res_cmp, res_cmp_arithmetic, data, value, num_ops, bitlength_data, bitlength_res, greater_than, equality);
}

void my::compare_with_eq_const_thread(int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[],
                                      int tid, uint8_t *res_cmp, uint8_t *res_eq, uint64_t *data, uint64_t value, int num_ops, int bitlength, bool greater_than)
{
    std::unique_ptr<MyFunctions> myfunction_ptr;
    if (tid & 1)
    {
        myfunction_ptr = std::make_unique<MyFunctions>(3 - party, iopackArr[tid], otpackArr[tid]);
    }
    else
    {
        myfunction_ptr = std::make_unique<MyFunctions>(party, iopackArr[tid], otpackArr[tid]);
    }
    myfunction_ptr->compare_with_eq_const(res_cmp, res_eq, data, value, num_ops, bitlength, greater_than);
}

void my::compare_with_eq_const_arithmetic_thread(int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[],
                                                 int tid, uint8_t *res_cmp, uint8_t *res_eq, uint64_t *res_cmp_arithmetic, uint64_t *res_eq_arithmetic, uint64_t *data, uint64_t value, int num_ops, int bitlength_data, int bitlength_res, bool greater_than)
{
    std::unique_ptr<MyFunctions> myfunction_ptr;
    if (tid & 1)
    {
        myfunction_ptr = std::make_unique<MyFunctions>(3 - party, iopackArr[tid], otpackArr[tid]);
    }
    else
    {
        myfunction_ptr = std::make_unique<MyFunctions>(party, iopackArr[tid], otpackArr[tid]);
    }
    myfunction_ptr->compare_with_eq_const_arithmetic(res_cmp, res_eq, res_cmp_arithmetic, res_eq_arithmetic, data, value, num_ops, bitlength_data, bitlength_res, greater_than);
}

void my::compare_thread(Compare_Type type, int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[], int tid,
                        uint64_t *data, uint64_t value, size_t num_ops, int bw_data, uint8_t *res_cmp, uint8_t *res_eq)
{
    std::unique_ptr<MyFunctions> myfunction_ptr;
    if (tid & 1)
    {
        myfunction_ptr = std::make_unique<MyFunctions>(3 - party, iopackArr[tid], otpackArr[tid]);
    }
    else
    {
        myfunction_ptr = std::make_unique<MyFunctions>(party, iopackArr[tid], otpackArr[tid]);
    }
    switch (type)
    {
    case LESS:
        myfunction_ptr->compare_const(res_cmp, data, value, num_ops, bw_data, false, false);
        break;
    case GREATER:
        myfunction_ptr->compare_const(res_cmp, data, value, num_ops, bw_data, true, false);
        break;
    case EQUAL:
        myfunction_ptr->check_equality(res_eq, data, value, num_ops, bw_data);
        break;
    case LESSwithEQUAL:
        myfunction_ptr->compare_with_eq_const(res_cmp, res_eq, data, value, num_ops, bw_data, false);
        break;
    case GREATERwithEQUAL:
        myfunction_ptr->compare_with_eq_const(res_cmp, res_eq, data, value, num_ops, bw_data, true);
        break;
    default:
        break;
    }
}

void my::compare_arithmetic_thread(Compare_Type type, int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[], int tid,
                                   uint64_t *data, uint64_t value, size_t num_ops, int bw_data, int bw_res,
                                   uint8_t *res_cmp, uint8_t *res_eq, uint64_t *res_cmp_arithmetic, uint64_t *res_eq_arithmetic)
{
    std::unique_ptr<MyFunctions> myfunction_ptr;
    if (tid & 1)
    {
        myfunction_ptr = std::make_unique<MyFunctions>(3 - party, iopackArr[tid], otpackArr[tid]);
    }
    else
    {
        myfunction_ptr = std::make_unique<MyFunctions>(party, iopackArr[tid], otpackArr[tid]);
    }
    switch (type)
    {
    case LESS:
        myfunction_ptr->compare_const_arithmetic(res_cmp, res_cmp_arithmetic, data, value, num_ops, bw_data, bw_res, false);
        break;
    case GREATER:
        myfunction_ptr->compare_const_arithmetic(res_cmp, res_cmp_arithmetic, data, value, num_ops, bw_data, bw_res, true);
        break;
    case EQUAL:
        myfunction_ptr->check_equality_arithmetic(res_eq, res_eq_arithmetic, data, value, num_ops, bw_data, bw_res);
        break;
    case LESSwithEQUAL:
        myfunction_ptr->compare_with_eq_const_arithmetic(res_cmp, res_eq, res_cmp_arithmetic, res_eq_arithmetic, data, value, num_ops, bw_data, bw_res, false);
        break;
    case GREATERwithEQUAL:
        myfunction_ptr->compare_with_eq_const_arithmetic(res_cmp, res_eq, res_cmp_arithmetic, res_eq_arithmetic, data, value, num_ops, bw_data, bw_res, true);
        break;
    default:
        break;
    }
};

void my::truncate_thread(Truncate_Type type, int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[], int tid,
                         int32_t dim, uint64_t *inA, uint64_t *outB, int32_t shift, int32_t bw, bool signed_arithmetic, uint8_t *msb_x)
{
    std::unique_ptr<Truncation> truncation_ptr;
    if (tid & 1)
    {
        truncation_ptr = std::make_unique<Truncation>(3 - party, iopackArr[tid], otpackArr[tid]);
    }
    else
    {
        truncation_ptr = std::make_unique<Truncation>(party, iopackArr[tid], otpackArr[tid]);
    }
    switch (type)
    {
    case TRUNCATE:
        truncation_ptr->truncate(dim, inA, outB, shift, bw, signed_arithmetic, msb_x);
        break;
    case DIV_POW2:
        truncation_ptr->div_pow2(dim, inA, outB, shift, bw, signed_arithmetic, msb_x);
        break;
    case TRUNCATE_RED_THEN_EXT:
        truncation_ptr->truncate_red_then_ext(dim, inA, outB, shift, bw, signed_arithmetic, msb_x);
        break;
    default:
        break;
    }
};

void my::mux_thread(int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[], int tid,
                    uint8_t *sel, uint64_t *x, uint64_t *y, int32_t size, int32_t bw_x, int32_t bw_y)
{
    std::unique_ptr<AuxProtocols> aux_ptr;
    if (tid & 1)
    {
        aux_ptr = std::make_unique<AuxProtocols>(3 - party, iopackArr[tid], otpackArr[tid]);
    }
    else
    {
        aux_ptr = std::make_unique<AuxProtocols>(party, iopackArr[tid], otpackArr[tid]);
    }
    aux_ptr->multiplexer(sel, x, y, size, bw_x, bw_y);
}

// signed_nm=false, compute_msnzb=true, 结果左移位s_dn+s_out-s_nm;
void my::secure_div_thread(int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[], int tid,
                           int32_t dim, uint64_t *nm, uint64_t *dn, uint64_t *out, int32_t bw_nm, int32_t bw_dn, int32_t bw_out,
                           int32_t s_nm, int32_t s_dn, int32_t s_out,
                           bool signed_nm, bool compute_msnzb)
{
    auto math_ptr = std::make_unique<MathFunctions>(party, iopackArr[tid], otpackArr[tid]);
    math_ptr->div(dim, nm, dn, out, bw_nm, bw_dn, bw_out, s_nm, s_dn, s_out, signed_nm, compute_msnzb);
}

void my::hadamard_product_thread(int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[], int tid,
                                 int32_t dim, uint64_t *inA, uint64_t *inB,
                                 uint64_t *outC, int32_t bwA, int32_t bwB,
                                 int32_t bwC, bool signed_arithmetic,
                                 bool signed_B, MultMode mode, uint8_t *msbA,
                                 uint8_t *msbB)
{
    auto linearOT_ptr = std::make_unique<LinearOT>(party, iopackArr[tid], otpackArr[tid]);
    linearOT_ptr->hadamard_product(dim, inA, inB, outC, bwA, bwB, bwC, signed_arithmetic, signed_B, mode, msbA, msbB);
}

void my::funcCompareWithEqConstArithmeticWrapper(int size, uint64_t *inp,
                                                 uint64_t *res_cmp_arithmetic, uint64_t *res_eq_arithmetic,
                                                 uint8_t *res_cmp, uint8_t *res_eq,
                                                 uint64_t compare_value, int bitlength_data, int bitlength_res, bool greater_than)
{
#ifdef MULTITHREADED_TRUNC
    std::thread truncThreads[num_threads];
    int chunk_size = size / num_threads;
    for (int i = 0; i < num_threads; i++)
    {
        int offset = i * chunk_size;
        int curSize;
        if (i == (num_threads - 1))
        {
            curSize = size - offset;
        }
        else
        {
            curSize = chunk_size;
        }
        // int curParty = party;
        // if (i & 1)
        //   curParty = 3 - curParty;
        truncThreads[i] = std::thread(
            my::compare_with_eq_const_arithmetic_thread,
            party, iopackArr, otpackArr, i, res_cmp + offset, res_eq + offset,
            res_cmp_arithmetic + offset,
            res_eq_arithmetic + offset,
            inp + offset, compare_value, curSize,
            bitlength_data, bitlength_res, greater_than);
    }
    for (int i = 0; i < num_threads; ++i)
    {
        truncThreads[i].join();
    }
#else
    std::cout << "" error !!!!!!!!!!! << std::endl;
//                           prg128Instance, size, inp, outp, divisor);
#endif
}

void my::ElemWiseVectorCompareWithEqArithmeticConst(
    int size, uint64_t *inp,
    uint64_t *res_cmp_arithmetic, uint64_t *res_eq_arithmetic,
    uint8_t *res_cmp, uint8_t *res_eq,
    uint64_t compare_value, int bitlength_data, int bitlength_res, bool greater_than)
{

#ifdef LOG_LAYERWISE
    INIT_TIMER;
    INIT_ALL_IO_DATA_SENT;
#endif

#ifdef SCI_OT
    my::funcCompareWithEqConstArithmeticWrapper(size, inp, res_cmp_arithmetic, res_eq_arithmetic, res_cmp, res_eq, compare_value, bitlength_data, bitlength_res, greater_than);
// #else
//   funcFieldDivWrapper(s1, inp, out, (uint64_t)divisor, nullptr);
#endif

#ifdef LOG_LAYERWISE
    auto temp = TIMER_TILL_NOW;
    CompareWithEqConstTimeInMilliSec += temp;
    std::cout << "Time in sec for current CompareWithEqConst = " << (temp / 1000.0)
              << std::endl;
    uint64_t curComm;
    FIND_ALL_IO_TILL_NOW(curComm);
    CompareWithEqConstCommSent += curComm;
#endif
    return;
}

void my::funcCompareWithEqConstWrapper(int size, uint64_t *inp,
                                       uint8_t *res_cmp, uint8_t *res_eq,
                                       uint64_t compare_value, int bitlength_data, bool greater_than)
{
#ifdef MULTITHREADED_TRUNC
    std::thread truncThreads[num_threads];
    int chunk_size = size / num_threads;
    for (int i = 0; i < num_threads; i++)
    {
        int offset = i * chunk_size;
        int curSize;
        if (i == (num_threads - 1))
        {
            curSize = size - offset;
        }
        else
        {
            curSize = chunk_size;
        }
        // int curParty = party;
        // if (i & 1)
        //   curParty = 3 - curParty;
        truncThreads[i] = std::thread(
            my::compare_with_eq_const_thread,
            party, iopackArr, otpackArr, i, res_cmp + offset, res_eq + offset,
            inp + offset, compare_value, curSize,
            bitlength_data, greater_than);
    }
    for (int i = 0; i < num_threads; ++i)
    {
        truncThreads[i].join();
    }
#else
    std::cout << "" error !!!!!!!!!!! << std::endl;
//                           prg128Instance, size, inp, outp, divisor);
#endif
}

void my::ElemWiseVectorCompareWithEqConst(
    int size, uint64_t *inp,
    uint8_t *res_cmp, uint8_t *res_eq,
    uint64_t compare_value, int bitlength_data, bool greater_than)
{

#ifdef LOG_LAYERWISE
    INIT_TIMER;
    INIT_ALL_IO_DATA_SENT;
#endif

#ifdef SCI_OT
    my::funcCompareWithEqConstWrapper(size, inp, res_cmp, res_eq, compare_value, bitlength_data, greater_than);
// #else
//   funcFieldDivWrapper(s1, inp, out, (uint64_t)divisor, nullptr);
#endif

#ifdef LOG_LAYERWISE
    auto temp = TIMER_TILL_NOW;
    CompareWithEqConstTimeInMilliSec += temp;
    std::cout << "Time in sec for current CompareWithEqConst = " << (temp / 1000.0)
              << std::endl;
    uint64_t curComm;
    FIND_ALL_IO_TILL_NOW(curComm);
    CompareWithEqConstCommSent += curComm;
#endif
    return;
}

void my::funcEqualArithmeticWrapper(int size, uint64_t *inp,
                                    uint64_t *res_eq_arithmetic,
                                    uint8_t *res_eq,
                                    uint64_t compare_value, int bitlength_data, int bitlength_res)
{
#ifdef MULTITHREADED_TRUNC
    std::thread truncThreads[num_threads];
    int chunk_size = size / num_threads;
    for (int i = 0; i < num_threads; i++)
    {
        int offset = i * chunk_size;
        int curSize;
        if (i == (num_threads - 1))
        {
            curSize = size - offset;
        }
        else
        {
            curSize = chunk_size;
        }
        // int curParty = party;
        // if (i & 1)
        //   curParty = 3 - curParty;
        truncThreads[i] = std::thread(
            my::equal_arithmetic_thread,
            party, iopackArr, otpackArr, i,
            res_eq + offset,
            res_eq_arithmetic + offset,
            inp + offset, compare_value, curSize,
            bitlength_data, bitlength_res);
    }
    for (int i = 0; i < num_threads; ++i)
    {
        truncThreads[i].join();
    }
#else
    std::cout << "" error !!!!!!!!!!! << std::endl;
//                           prg128Instance, size, inp, outp, divisor);
#endif
}

void my::ElemWiseVectorEqualArithmetic(int size, uint64_t *inp,
                                       uint64_t *res_eq_arithmetic,
                                       uint8_t *res_eq,
                                       uint64_t compare_value, int bitlength_data, int bitlength_res)
{

#ifdef LOG_LAYERWISE
    INIT_TIMER;
    INIT_ALL_IO_DATA_SENT;
#endif

#ifdef SCI_OT
    my::funcEqualArithmeticWrapper(size, inp, res_eq_arithmetic, res_eq, compare_value, bitlength_data, bitlength_res);
// #else
//   funcFieldDivWrapper(s1, inp, out, (uint64_t)divisor, nullptr);
#endif

#ifdef LOG_LAYERWISE
    auto temp = TIMER_TILL_NOW;
    EqualTimeInMilliSec += temp;
    std::cout << "Time in sec for current Equal = " << (temp / 1000.0)
              << std::endl;
    uint64_t curComm;
    FIND_ALL_IO_TILL_NOW(curComm);
    EqualCommSent += curComm;
#endif
    return;
}

void my::funcEqualWrapper(int size, uint64_t *inp,
                          uint8_t *res_eq,
                          uint64_t compare_value,
                          int bitlength_data)
{
#ifdef MULTITHREADED_TRUNC
    std::thread truncThreads[num_threads];
    int chunk_size = size / num_threads;
    for (int i = 0; i < num_threads; i++)
    {
        int offset = i * chunk_size;
        int curSize;
        if (i == (num_threads - 1))
        {
            curSize = size - offset;
        }
        else
        {
            curSize = chunk_size;
        }
        // int curParty = party;
        // if (i & 1)
        //   curParty = 3 - curParty;
        truncThreads[i] = std::thread(
            my::equal_thread,
            party, iopackArr, otpackArr, i,
            res_eq + offset,
            inp + offset, compare_value, curSize,
            bitlength_data);
    }
    for (int i = 0; i < num_threads; ++i)
    {
        truncThreads[i].join();
    }
#else
    std::cout << "" error !!!!!!!!!!! << std::endl;
//                           prg128Instance, size, inp, outp, divisor);
#endif
}

void my::ElemWiseVectorEqual(int size, uint64_t *inp,
                             uint8_t *res_eq,
                             uint64_t compare_value,
                             int bitlength_data)
{

#ifdef LOG_LAYERWISE
    INIT_TIMER;
    INIT_ALL_IO_DATA_SENT;
#endif

#ifdef SCI_OT
    my::funcEqualWrapper(size, inp, res_eq, compare_value, bitlength_data);
// #else
//   funcFieldDivWrapper(s1, inp, out, (uint64_t)divisor, nullptr);
#endif

#ifdef LOG_LAYERWISE
    auto temp = TIMER_TILL_NOW;
    EqualTimeInMilliSec += temp;
    std::cout << "Time in sec for current Equal = " << (temp / 1000.0)
              << std::endl;
    uint64_t curComm;
    FIND_ALL_IO_TILL_NOW(curComm);
    EqualCommSent += curComm;
#endif
    return;
}

void my::funcCompareConstArithmeticWrapper(int size,
                                           uint64_t *inp,
                                           uint64_t *res_cmp_arithmetic,
                                           uint8_t *res_cmp,
                                           uint64_t compare_value,
                                           int bitlength_data,
                                           int bitlength_res,
                                           bool greater_than,
                                           bool equality)
{
#ifdef MULTITHREADED_TRUNC
    std::thread truncThreads[num_threads];
    int chunk_size = size / num_threads;
    for (int i = 0; i < num_threads; i++)
    {
        int offset = i * chunk_size;
        int curSize;
        if (i == (num_threads - 1))
        {
            curSize = size - offset;
        }
        else
        {
            curSize = chunk_size;
        }
        // int curParty = party;
        // if (i & 1)
        //   curParty = 3 - curParty;
        truncThreads[i] = std::thread(
            my::compare_const_arithmetic_thread,
            party, iopackArr, otpackArr, i, res_cmp + offset,
            res_cmp_arithmetic + offset,
            inp + offset, compare_value, curSize,
            bitlength_data, bitlength_res, greater_than, equality);
    }
    for (int i = 0; i < num_threads; ++i)
    {
        truncThreads[i].join();
    }
#else
    std::cout << "" error !!!!!!!!!!! << std::endl;
//                           prg128Instance, size, inp, outp, divisor);
#endif
}

void my::ElemWiseVectorCompareConstArithmetic(int size,
                                              uint64_t *inp,
                                              uint64_t *res_cmp_arithmetic,
                                              uint8_t *res_cmp,
                                              uint64_t compare_value,
                                              int bitlength_data,
                                              int bitlength_res,
                                              bool greater_than)
{

#ifdef LOG_LAYERWISE
    INIT_TIMER;
    INIT_ALL_IO_DATA_SENT;
#endif

#ifdef SCI_OT
    my::funcCompareConstArithmeticWrapper(size, inp, res_cmp_arithmetic, res_cmp, compare_value, bitlength_data, bitlength_res, greater_than);
#endif

#ifdef LOG_LAYERWISE
    auto temp = TIMER_TILL_NOW;
    CompareConstTimeInMilliSec += temp;
    std::cout << "Time in sec for current CompareConst = " << (temp / 1000.0)
              << std::endl;
    uint64_t curComm;
    FIND_ALL_IO_TILL_NOW(curComm);
    CompareConstCommSent += curComm;
#endif
    return;
}

void my::funcCompareConstWrapper(int size,
                                 uint64_t *inp,
                                 uint8_t *res_cmp,
                                 uint64_t compare_value,
                                 int bitlength_data,
                                 bool greater_than,
                                 bool equality)
{
#ifdef MULTITHREADED_TRUNC
    std::thread truncThreads[num_threads];
    int chunk_size = size / num_threads;
    for (int i = 0; i < num_threads; i++)
    {
        int offset = i * chunk_size;
        int curSize;
        if (i == (num_threads - 1))
        {
            curSize = size - offset;
        }
        else
        {
            curSize = chunk_size;
        }
        // int curParty = party;
        // if (i & 1)
        //   curParty = 3 - curParty;
        truncThreads[i] = std::thread(
            my::compare_const_thread,
            party, iopackArr, otpackArr, i, res_cmp + offset,
            inp + offset, compare_value, curSize,
            bitlength_data, greater_than, equality);
    }
    for (int i = 0; i < num_threads; ++i)
    {
        truncThreads[i].join();
    }
#else
    std::cout << "" error !!!!!!!!!!! << std::endl;
//                           prg128Instance, size, inp, outp, divisor);
#endif
}

void my::ElemWiseVectorCompareConst(int size,
                                    uint64_t *inp,
                                    uint8_t *res_cmp,
                                    uint64_t compare_value,
                                    int bitlength_data,
                                    bool greater_than)
{

#ifdef LOG_LAYERWISE
    INIT_TIMER;
    INIT_ALL_IO_DATA_SENT;
#endif

#ifdef SCI_OT
    my::funcCompareConstWrapper(size, inp, res_cmp, compare_value, bitlength_data, greater_than);
#endif

#ifdef LOG_LAYERWISE
    auto temp = TIMER_TILL_NOW;
    CompareConstTimeInMilliSec += temp;
    std::cout << "Time in sec for current CompareConst = " << (temp / 1000.0)
              << std::endl;
    uint64_t curComm;
    FIND_ALL_IO_TILL_NOW(curComm);
    CompareConstCommSent += curComm;
#endif
    return;
}

void my::SecureDiv_Thread(
    int32_t tid,
    int32_t dim,
    // numerator
    uint64_t *nm,
    // denominator
    uint64_t *dn,
    // output
    uint64_t *out,
    // bitwidths
    int32_t bw_nm, int32_t bw_dn, int32_t bw_out,
    // scales
    int32_t s_nm, int32_t s_dn, int32_t s_out,
    bool signed_nm,
    bool compute_msnzb)
{
    mathArr[tid]->div(dim, nm, dn, out, bw_nm, bw_dn, bw_out, s_nm, s_dn, s_out, signed_nm, compute_msnzb);
}

void my::funcSecureDivWrapper(
    int size,
    // numerator
    uint64_t *nm,
    // denominator
    uint64_t *dn,
    // output
    uint64_t *out,
    // bitwidths
    int32_t bw_nm, int32_t bw_dn, int32_t bw_out,
    // scales
    int32_t s_nm, int32_t s_dn, int32_t s_out,
    bool signed_nm,
    bool compute_msnzb)
{
#ifdef MULTITHREADED_TRUNC
    std::thread truncThreads[num_threads];
    int chunk_size = size / num_threads;
    for (int i = 0; i < num_threads; i++)
    {
        int offset = i * chunk_size;
        int curSize;
        if (i == (num_threads - 1))
        {
            curSize = size - offset;
        }
        else
        {
            curSize = chunk_size;
        }
        truncThreads[i] = std::thread(
            my::SecureDiv_Thread,
            i, curSize, nm + offset, dn + offset, out + offset, bw_nm, bw_dn, bw_out, s_nm, s_dn, s_out, signed_nm, compute_msnzb);
    }
    for (int i = 0; i < num_threads; ++i)
    {
        truncThreads[i].join();
    }
#else
    std::cout << "" error !!!!!!!!!!! << std::endl;
//                           prg128Instance, size, inp, outp, divisor);
#endif
}

void my::ElemWiseVectorSecureDiv(
    int size,
    // numerator
    uint64_t *nm,
    // denominator
    uint64_t *dn,
    // output
    uint64_t *out,
    // bitwidths
    int32_t bw_nm, int32_t bw_dn, int32_t bw_out,
    // scales
    int32_t s_nm, int32_t s_dn, int32_t s_out,
    bool signed_nm,
    bool compute_msnzb)
{

#ifdef LOG_LAYERWISE
    INIT_TIMER;
    INIT_ALL_IO_DATA_SENT;
#endif

#ifdef SCI_OT
    my::funcSecureDivWrapper(size, nm, dn, out, bw_nm, bw_dn, bw_out, s_nm, s_dn, s_out, signed_nm, compute_msnzb);
#endif

#ifdef LOG_LAYERWISE
    auto temp = TIMER_TILL_NOW;
    SecureDivTimeInMilliSec += temp;
    std::cout << "Time in sec for current SecureDiv = " << (temp / 1000.0)
              << std::endl;
    uint64_t curComm;
    FIND_ALL_IO_TILL_NOW(curComm);
    SecureDivCommSent += curComm;
#endif
    return;
}

void my::VectorSortEmulate(int size, uint64_t *inp, uint64_t *res, int bitlength)
{

#ifdef LOG_LAYERWISE
    INIT_TIMER;
    INIT_ALL_IO_DATA_SENT;
#endif

#ifdef SCI_OT
    int tmp = ceil(log2(size));
    int round = tmp * (tmp + 1) / 2;
    size_t compare_size = size / 2;
    auto mask = my::mask_gen_uint64_t(bitlength);

    vector<uint64_t> input_vector(inp, inp + size);

    vector<uint64_t> cmp_input_1(compare_size);
    vector<uint64_t> cmp_input_2(compare_size);
    vector<uint64_t> cmp_input(compare_size);
    vector<uint8_t> cmp_res_uint8(compare_size);

    vector<uint64_t> mux_result(compare_size);
    vector<uint64_t> cmp_reslut(compare_size);

    uint64_t compare_value = 0;
    for (int i = 0; i < round; i++)
    {
        for(size_t i = 0; i<compare_size; i++){
            cmp_input_1[i] = input_vector[2*i];
            cmp_input_2[i] = input_vector[2*i+1];
            cmp_input[i] = (cmp_input_1[i]-cmp_input_2[i])&mask;
        }
        my::funcCompareConstWrapper(compare_size, cmp_input.data(), cmp_res_uint8.data(), compare_value, bitlength);
        my::funcMUXWrapper(cmp_res_uint8.data(),cmp_input.data(),mux_result.data(),compare_size,bitlength,bitlength);
        
        for(size_t i = 0; i<compare_size; i++){
            cmp_reslut[i] = (mux_result[i] + cmp_input_2[i])& mask;
            input_vector[2*i] = cmp_reslut[i];
            input_vector[2*i+1] = (cmp_input_1[i] + cmp_input_2[i]) - cmp_reslut[i];
        }
    }
#endif

#ifdef LOG_LAYERWISE
    auto temp = TIMER_TILL_NOW;
    SortEmulateTimeInMilliSec += temp;
    std::cout << "Time in sec for current SortEmulate = " << (temp / 1000.0) << " compare times = " << round << std::endl;
    uint64_t curComm;
    FIND_ALL_IO_TILL_NOW(curComm);
    SortEmulateCommSent += curComm;
#endif
    return;
}

void my::MUX_Thread(int32_t tid, uint8_t *sel, uint64_t *x, uint64_t *y,
                    int32_t size, int32_t bw_x, int32_t bw_y)
{
    auxArr[tid]->multiplexer(sel, x, y, size, bw_x, bw_y);
}

void my::funcMUXWrapper(uint8_t *sel, uint64_t *x, uint64_t *y,
                        int32_t size, int32_t bw_x, int32_t bw_y)
{
#ifdef MULTITHREADED_TRUNC
    std::thread truncThreads[num_threads];
    int chunk_size = size / num_threads;
    for (int i = 0; i < num_threads; i++)
    {
        int offset = i * chunk_size;
        int curSize;
        if (i == (num_threads - 1))
        {
            curSize = size - offset;
        }
        else
        {
            curSize = chunk_size;
        }
        truncThreads[i] = std::thread(
            my::MUX_Thread, i, sel + offset, x + offset, y + offset, curSize, bw_x, bw_y);
    }
    for (int i = 0; i < num_threads; ++i)
    {
        truncThreads[i].join();
    }
#else
    std::cout << "" error !!!!!!!!!!! << std::endl;
//                           prg128Instance, size, inp, outp, divisor);
#endif
}

void my::ElemWiseVectorMUX(uint8_t *sel, uint64_t *x, uint64_t *y,
                           int32_t size, int32_t bw_x, int32_t bw_y)
{

#ifdef LOG_LAYERWISE
    INIT_TIMER;
    INIT_ALL_IO_DATA_SENT;
#endif

#ifdef SCI_OT
    my::funcMUXWrapper(sel, x, y, size, bw_x, bw_y);
#endif

#ifdef LOG_LAYERWISE
    auto temp = TIMER_TILL_NOW;
    MUXTimeInMilliSec += temp;
    std::cout << "Time in sec for current MUX = " << (temp / 1000.0)
              << std::endl;
    uint64_t curComm;
    FIND_ALL_IO_TILL_NOW(curComm);
    MUXCommSent += curComm;
#endif
    return;
}

void my::ElemWiseVectorMUXRounds(int32_t round, uint8_t *sel, uint64_t *x, uint64_t *y,
                                 int32_t size, int32_t bw_x, int32_t bw_y)
{

#ifdef LOG_LAYERWISE
    INIT_TIMER;
    INIT_ALL_IO_DATA_SENT;
#endif

#ifdef SCI_OT
    for (int i = 0; i < round; i++)
    {
        my::funcMUXWrapper(sel, x, y, size, bw_x, bw_y);
    }
#endif

#ifdef LOG_LAYERWISE
    auto temp = TIMER_TILL_NOW;
    MUXRoundsTimeInMilliSec += temp;
    std::cout << "Time in sec for current MUXRounds = " << (temp / 1000.0)
              << std::endl;
    uint64_t curComm;
    FIND_ALL_IO_TILL_NOW(curComm);
    MUXRoundsCommSent += curComm;
#endif
    return;
}

void my::MultipleLUT_Thread(int32_t tid, uint64_t **spec, uint64_t *x, uint64_t *y,
                            int32_t size, int32_t bw_x, int32_t bw_y)
{
    auxArr[tid]->lookup_table(spec, x, y, size, bw_x, bw_y);
}

void my::funcMultipleLUTWrapper(uint64_t **spec, uint64_t *x, uint64_t *y,
                                int32_t size, int32_t bw_x, int32_t bw_y)
{
#ifdef MULTITHREADED_TRUNC
    std::thread truncThreads[num_threads];
    int chunk_size = size / num_threads;
    for (int i = 0; i < num_threads; i++)
    {
        int offset = i * chunk_size;
        int curSize;
        if (i == (num_threads - 1))
        {
            curSize = size - offset;
        }
        else
        {
            curSize = chunk_size;
        }
        truncThreads[i] = std::thread(
            my::MultipleLUT_Thread, i, spec + offset, x + offset, y + offset, curSize, bw_x, bw_y);
    }
    for (int i = 0; i < num_threads; ++i)
    {
        truncThreads[i].join();
    }
#else
    std::cout << "" error !!!!!!!!!!! << std::endl;
//                           prg128Instance, size, inp, outp, divisor);
#endif
}

void my::MultipleLUT(uint64_t **spec, uint64_t *x, uint64_t *y,
                     int32_t size, int32_t bw_x, int32_t bw_y)
{

#ifdef LOG_LAYERWISE
    INIT_TIMER;
    INIT_ALL_IO_DATA_SENT;
#endif

#ifdef SCI_OT
    my::funcMultipleLUTWrapper(spec, x, y, size, bw_x, bw_y);
#endif

#ifdef LOG_LAYERWISE
    auto temp = TIMER_TILL_NOW;
    MultipleLUTTimeInMilliSec += temp;
    std::cout << "Time in sec for current MultipleLUT = " << (temp / 1000.0)
              << std::endl;
    uint64_t curComm;
    FIND_ALL_IO_TILL_NOW(curComm);
    MultipleLUTCommSent += curComm;
#endif
    return;
}
