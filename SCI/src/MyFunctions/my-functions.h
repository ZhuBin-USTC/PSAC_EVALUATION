#ifndef _MYFUNCTION_H
#define _MYFUNCTION_H

#include "BuildingBlocks/aux-protocols.h"
#include "Millionaire/millionaire.h"
#include "Millionaire/equality.h"
#include "library_fixed.h"
#include <execution>

#define MULTITHREADED_TRUNC

using ShareA = uint64_t;
using ShareB = uint8_t;

enum Truncate_Type
{
    TRUNCATE,
    DIV_POW2,
    TRUNCATE_RED_THEN_EXT
};

enum Compare_Type
{
    LESS,
    GREATER,
    EQUAL,
    LESSwithEQUAL,
    GREATERwithEQUAL
};

class MyFunctions
{
public:
    int party;
    sci::IOPack *iopack;
    sci::OTPack *otpack;
    // Equality *equality;
    // AuxProtocols *aux;
    std::unique_ptr<AuxProtocols> aux;
    std::unique_ptr<Equality> equality;
    std::unique_ptr<MillionaireProtocol> millionaire;

    MyFunctions(int party, sci::IOPack *iopack, sci::OTPack *otpack);
    // ~MyFunctions();

    void check_equality(uint8_t *res_eq, uint64_t *data, uint64_t value, int num_eqs, int bitlength, int radix_base = MILL_PARAM);
    void check_equality_arithmetic(uint8_t *res_eq, uint64_t *res_eq_arithmetic, uint64_t *data, uint64_t value, int num_eqs, int bitlength_data, int bitlength_res_arithmetic, int radix_base = MILL_PARAM);

    void compare_const(uint8_t *res_cmp, uint64_t *data, uint64_t value, int num_cmps, int bitlength_data, bool greater_than = true, bool equality = false, int radix_base = MILL_PARAM);
    void compare_const_arithmetic(uint8_t *res_cmp, uint64_t *res_cmp_arithmetic, uint64_t *data, uint64_t value, int num_cmps, int bitlength_data, int bitlength_res_arithmetic, bool greater_than = true, bool equality = false, int radix_base = MILL_PARAM);

    // void MSB_with_eq(uint8_t *res_cmp, uint8_t *res_eq, uint64_t *data, uint64_t value, int num_cmps, int bitlength, bool greater_than = true, int radix_base = MILL_PARAM);
    void compare_with_eq_const(uint8_t *res_cmp, uint8_t *res_eq, uint64_t *data, uint64_t value, int num_cmps, int bitlength_data, bool greater_than = true, int radix_base = MILL_PARAM);
    void compare_with_eq_const_arithmetic(uint8_t *res_cmp, uint8_t *res_eq, uint64_t *res_cmp_arithmetic, uint64_t *res_eq_arithmetic, uint64_t *data, uint64_t value, int num_cmps, int bitlength_data, int bitlength_res_arithmetic, bool greater_than = true, int radix_base = MILL_PARAM);

    void softmax(int32_t dim, uint64_t *x, uint64_t *y, uint64_t max_x, int32_t bw_x,
                 int32_t bw_y, int32_t s_x, int32_t s_y);
};

namespace my
{
    ShareA shareA_const(int party, uint64_t x);
    ShareB shareB_const(int party, uint8_t x);

    void set_value_int(int &var, int value);
    void set_value_uint64_t(uint64_t &var, uint64_t value);
    void set_value_double(double &var, double value);

    void setup_io_and_ots(int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[], int num_threads, int port, std::string address);

    std::vector<ShareA> append_vectors(const std::vector<ShareA> input1, const std::vector<ShareA> input2);

    uint64_t mask_gen_uint64_t(int bw);

    std::vector<uint64_t> reconstruct_uint64_t(int party, int reconstruct_party, sci::IOPack *iopackArr[], std::vector<uint64_t> &sharing, int bw);
    std::vector<uint8_t> reconstruct_uint8_t(int party, int reconstruct_party, sci::IOPack *iopackArr[], std::vector<uint8_t> &sharing);

    // 专门用于本地测试的，不加HE协议
    std::unique_ptr<uint64_t[]> generate_and_send_sharing_vector_plaintext(sci::IOPack *iopack, std::unique_ptr<uint64_t[]> &data_vector_ptr, int dim, int bitlength);
    std::unique_ptr<uint64_t[]> recv_sharing_vector_plaintext(sci::IOPack *iopack, int dim);
    std::unique_ptr<uint64_t[]> generated_worker_data_arithmetic(int num_of_data, int bitlength_arithmetic, uint64_t min, uint64_t max);
    std::unique_ptr<uint64_t[]> data_arithmetic_to_frequency(int dim_frequency, uint64_t left, uint64_t step, std::unique_ptr<uint64_t[]> &data_arithmetic_ptr, int num_of_data);

    std::vector<ShareA> generate_and_send_sharing_vector_plaintext(sci::IOPack *iopack, std::vector<uint64_t> &data_vector, int bitlength);
    std::vector<uint64_t> recv_sharing_vector_plaintext_vector(sci::IOPack *iopack, int dim);
    std::vector<uint64_t> generate_numeric_data_vector(int num_of_data, int bitlength_arithmetic, uint64_t min, uint64_t max);
    std::vector<uint64_t> numeric_data_to_hist(int dim_frequency, uint64_t left, uint64_t step, std::vector<uint64_t> &data_arithmetic, int num_of_data);

    void equal_thread(int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[],
                      int tid, uint8_t *res_eq, uint64_t *data, uint64_t value, int num_ops, int bitlength);

    void equal_arithmetic_thread(int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[],
                                 int tid, uint8_t *res_eq, uint64_t *res_eq_arithmetic, uint64_t *data, uint64_t value, int num_ops, int bitlength_data, int bitlength_res);

    void compare_const_thread(int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[],
                              int tid, uint8_t *res_cmp, uint64_t *data, uint64_t value, int num_ops, int bitlength, bool greater_than = false, bool equality = false);

    void compare_const_arithmetic_thread(int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[],
                                         int tid, uint8_t *res_cmp, uint64_t *res_cmp_arithmetic, uint64_t *data, uint64_t value, int num_ops, int bitlength_data, int bitlength_res, bool greater_than = false, bool equality = false);

    void compare_with_eq_const_thread(int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[],
                                      int tid, uint8_t *res_cmp, uint8_t *res_eq, uint64_t *data, uint64_t value, int num_ops, int bitlength, bool greater_than = false);

    void compare_with_eq_const_arithmetic_thread(int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[],
                                                 int tid, uint8_t *res_cmp, uint8_t *res_eq, uint64_t *res_cmp_arithmetic, uint64_t *res_eq_arithmetic, uint64_t *data, uint64_t value, int num_ops, int bitlength_data, int bitlength_res, bool greater_than = false);

    void compare_thread(Compare_Type type, int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[], int tid,
                        uint64_t *data, uint64_t value, size_t num_ops, int bw_data, uint8_t *res_cmp = nullptr, uint8_t *res_eq = nullptr);

    void compare_arithmetic_thread(Compare_Type type, int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[], int tid,
                                   uint64_t *data, uint64_t value, size_t num_ops, int bw_data, int bw_res, uint8_t *res_cmp = nullptr, uint8_t *res_eq = nullptr, uint64_t *res_cmp_arithmetic = nullptr, uint64_t *res_eq_arithmetic = nullptr);

    void truncate_thread(Truncate_Type type, int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[], int tid,
                         // Size of vector
                         int32_t dim,
                         // input vector
                         uint64_t *inA,
                         // output vector
                         uint64_t *outB,
                         // right shift amount
                         int32_t shift,
                         // Input and output bitwidth
                         int32_t bw,
                         // signed truncation?
                         bool signed_arithmetic = true,
                         // msb of input vector elements
                         uint8_t *msb_x = nullptr);

    void mux_thread(int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[], int tid,
                    uint8_t *sel, uint64_t *x, uint64_t *y, int32_t size, int32_t bw_x, int32_t bw_y);

    void secure_div_thread(int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[], int tid,
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
                           bool compute_msnzb);

    void hadamard_product_thread(int party, sci::IOPack *iopackArr[], sci::OTPack *otpackArr[], int tid,
                                 int32_t dim,
                                 // input share vector
                                 uint64_t *inA, uint64_t *inB,
                                 // output share vector
                                 uint64_t *outC,
                                 // bitwidths
                                 int32_t bwA, int32_t bwB, int32_t bwC,
                                 bool signed_arithmetic = true,
                                 // take B as signed input?
                                 bool signed_B = true, MultMode mode = MultMode::None,
                                 uint8_t *msbA = nullptr, uint8_t *msbB = nullptr);

    // todo zhubin 08
    void funcCompareWithEqConstArithmeticWrapper(
        int size,
        uint64_t *inp,
        uint64_t *res_cmp_arithmetic,
        uint64_t *res_eq_arithmetic,
        uint8_t *res_cmp,
        uint8_t *res_eq,
        uint64_t compare_value,
        int bitlength_data,
        int bitlength_res,
        bool greater_than = true);
    void ElemWiseVectorCompareWithEqArithmeticConst(
        int size,
        uint64_t *inp,
        uint64_t *res_cmp_arithmetic,
        uint64_t *res_eq_arithmetic,
        uint8_t *res_cmp,
        uint8_t *res_eq,
        uint64_t compare_value,
        int bitlength_data,
        int bitlength_res,
        bool greater_than = true);

    void funcCompareWithEqConstWrapper(
        int size, uint64_t *inp,
        uint8_t *res_cmp, uint8_t *res_eq,
        uint64_t compare_value, int bitlength_data, bool greater_than);

    void ElemWiseVectorCompareWithEqConst(
        int size, uint64_t *inp,
        uint8_t *res_cmp, uint8_t *res_eq,
        uint64_t compare_value, int bitlength_data, bool greater_than);

    void funcEqualArithmeticWrapper(
        int size,
        uint64_t *inp,
        uint64_t *res_eq_arithmetic,
        uint8_t *res_eq,
        uint64_t compare_value,
        int bitlength_data,
        int bitlength_res);
    void ElemWiseVectorEqualArithmetic(
        int size,
        uint64_t *inp,
        uint64_t *res_eq_arithmetic,
        uint8_t *res_eq,
        uint64_t compare_value,
        int bitlength_data,
        int bitlength_res);

    void funcEqualWrapper(
        int size,
        uint64_t *inp,
        uint8_t *res_eq,
        uint64_t compare_value,
        int bitlength_data);
    void ElemWiseVectorEqual(
        int size,
        uint64_t *inp,
        uint8_t *res_eq,
        uint64_t compare_value,
        int bitlength_data);

    void funcCompareConstArithmeticWrapper(
        int size,
        uint64_t *inp,
        uint64_t *res_cmp_arithmetic,
        uint8_t *res_cmp,
        uint64_t compare_value,
        int bitlength_data,
        int bitlength_res,
        bool greater_than = false,
        bool equality = false);
    void ElemWiseVectorCompareConstArithmetic(
        int size,
        uint64_t *inp,
        uint64_t *res_cmp_arithmetic,
        uint8_t *res_cmp,
        uint64_t compare_value,
        int bitlength_data,
        int bitlength_res,
        bool greater_than = false);

    void funcCompareConstWrapper(
        int size,
        uint64_t *inp,
        uint8_t *res_cmp,
        uint64_t compare_value,
        int bitlength_data,
        bool greater_than = false,
        bool equality = false);
    void ElemWiseVectorCompareConst(
        int size,
        uint64_t *inp,
        uint8_t *res_cmp,
        uint64_t compare_value,
        int bitlength_data,
        bool greater_than = false);

    void SecureDiv_Thread(
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
        bool signed_nm = false,
        bool compute_msnzb = true);
    void funcSecureDivWrapper(
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
        bool signed_nm = false,
        bool compute_msnzb = true);
    void ElemWiseVectorSecureDiv(
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
        bool signed_nm = false,
        bool compute_msnzb = true);

    void VectorSortEmulate(int size, uint64_t *inp, uint64_t *res, int bitlength);

    void MUX_Thread(int32_t tid, uint8_t *sel, uint64_t *x, uint64_t *y,
                    int32_t size, int32_t bw_x, int32_t bw_y);

    void funcMUXWrapper(uint8_t *sel, uint64_t *x, uint64_t *y,
                        int32_t size, int32_t bw_x, int32_t bw_y);
    void ElemWiseVectorMUX(uint8_t *sel, uint64_t *x, uint64_t *y,
                           int32_t size, int32_t bw_x, int32_t bw_y);
    void ElemWiseVectorMUXRounds(int32_t round, uint8_t *sel, uint64_t *x, uint64_t *y,
                                 int32_t size, int32_t bw_x, int32_t bw_y);

    void MultipleLUT_Thread(int32_t tid, uint64_t **spec, uint64_t *x, uint64_t *y,
                            int32_t size, int32_t bw_x, int32_t bw_y);
    void funcMultipleLUTWrapper(uint64_t **spec, uint64_t *x, uint64_t *y,
                                int32_t size, int32_t bw_x, int32_t bw_y);
    void MultipleLUT(uint64_t **spec, uint64_t *x, uint64_t *y,
                     int32_t size, int32_t bw_x, int32_t bw_y);
}

#endif