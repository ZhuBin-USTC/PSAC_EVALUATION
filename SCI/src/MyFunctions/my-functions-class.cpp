#include "MyFunctions/my-functions.h"

MyFunctions::MyFunctions(int party, sci::IOPack *iopack, sci::OTPack *otpack)
{
    this->party = party;
    this->iopack = iopack;
    this->otpack = otpack;
    this->equality = std::make_unique<Equality>(party, iopack, otpack);
    this->aux = std::make_unique<AuxProtocols>(party, iopack, otpack);
    this->millionaire = std::make_unique<MillionaireProtocol>(party, iopack, otpack);
}
// MyFunctions::~MyFunctions()
// {
// }

// 比较大小，为0表示yes，为1表示no

void MyFunctions::check_equality(uint8_t *res_eq, uint64_t *data, uint64_t value, int num_eqs, int bitlength, int radix_base)
{
    uint64_t mask_data = (bitlength == 64 ? -1 : ((1ULL << bitlength) - 1));
    uint64_t value_masked = value & mask_data;
    if (party == sci::ALICE)
    {
        // uint64_t *data_tmp = new uint64_t[num_eqs];
        auto data_tmp_ptr = std::make_unique<uint64_t[]>(num_eqs);
        for (int i = 0; i < num_eqs; i++)
        {
            data_tmp_ptr[i] = (mask_data + 1 - data[i] + value_masked) & mask_data;
        }
        equality->check_equality(res_eq, data_tmp_ptr.get(), num_eqs, bitlength, radix_base);
    }
    else if (party == sci::BOB)
    {
        equality->check_equality(res_eq, data, num_eqs, bitlength, radix_base);
    }
}

// void MyFunctions::compare_const(uint8_t *res_cmp, uint64_t *data, uint64_t value, int num_cmps, int bitlength, bool greater_than, bool equality, int radix_base)
// {
//     uint64_t mask_data = (bitlength == 64 ? -1 : ((1ULL << bitlength) - 1));
//     uint64_t value_masked = value & mask_data;
//     if (party == sci::ALICE)
//     {
//         auto data_tmp_ptr = std::make_unique<uint64_t[]>(num_cmps);
//         for (int i = 0; i < num_cmps; i++)
//         {
//             data_tmp_ptr[i] = (mask_data + 1 - data[i] + value_masked) & mask_data;
//         }
//         millionaire->compare(res_cmp, data_tmp_ptr.get(), num_cmps, bitlength, greater_than, equality, radix_base);
//     }
//     else if (party == sci::BOB)
//     {
//         millionaire->compare(res_cmp, data, num_cmps, bitlength, greater_than, equality, radix_base);
//     }
// }

void MyFunctions::compare_const(uint8_t *res_cmp, uint64_t *data, uint64_t value, int num_cmps, int bitlength_data, bool greater_than, bool equality, int radix_base)
{
    assert(bitlength_data <= 64);

    int32_t shift = bitlength_data - 1;
    uint64_t shift_mask = (shift == 64 ? -1 : ((1ULL << shift) - 1));
    uint64_t mask_data = (bitlength_data == 64 ? -1 : ((1ULL << bitlength_data) - 1));

    auto tmp_x = std::make_unique<uint64_t[]>(num_cmps);
    auto msb_xb = std::make_unique<uint8_t[]>(num_cmps);

    if (party == sci::ALICE)
    {
        for (int i = 0; i < num_cmps; i++)
        {
            // 原始数据减去value+1才行，然后读取MSB符号位。这里由Alice做
            uint64_t data_tmp = (data[i] - (value)) & mask_data;
            // uint64_t data_tmp = (data[i] - (value + 1)) & mask_data;
            tmp_x[i] = data_tmp & shift_mask;
            msb_xb[i] = (data_tmp >> shift) & 1;
        }
    }
    else if (party == sci::BOB)
    {
        for (int i = 0; i < num_cmps; i++)
        {
            tmp_x[i] = data[i] & shift_mask;
            msb_xb[i] = (data[i] >> shift) & 1;
            tmp_x[i] = (shift_mask - tmp_x[i]) & shift_mask;
        }
    }
    millionaire->compare(res_cmp, tmp_x.get(), num_cmps, bitlength_data - 1, true, equality, radix_base);

    for (int i = 0; i < num_cmps; i++)
    {
        res_cmp[i] = res_cmp[i] ^ msb_xb[i];
    }
    if (greater_than == true && party == sci::ALICE)
    {
        for (int i = 0; i < num_cmps; i++)
        {
            res_cmp[i] = res_cmp[i] ^ 1;
        }
    }
}

// less than
// 从MSB来的
void MyFunctions::compare_with_eq_const(uint8_t *res_cmp, uint8_t *res_eq, uint64_t *data, uint64_t value, int num_cmps, int bitlength_data, bool greater_than, int radix_base)
{
    assert(bitlength_data <= 64);

    int32_t shift = bitlength_data - 1;
    uint64_t shift_mask = (shift == 64 ? -1 : ((1ULL << shift) - 1));
    uint64_t mask_data = (bitlength_data == 64 ? -1 : ((1ULL << bitlength_data) - 1));

    auto tmp_x = std::make_unique<uint64_t[]>(num_cmps);
    auto msb_xb = std::make_unique<uint8_t[]>(num_cmps);

    if (party == sci::ALICE)
    {
        for (int i = 0; i < num_cmps; i++)
        {
            // 原始数据减去value+1才行，然后读取MSB符号位。这里由Alice做
            uint64_t data_tmp = (data[i] - (value + 1)) & mask_data;
            tmp_x[i] = data_tmp & shift_mask;
            msb_xb[i] = (data_tmp >> shift) & 1;
        }
    }
    else if (party == sci::BOB)
    {
        for (int i = 0; i < num_cmps; i++)
        {
            tmp_x[i] = data[i] & shift_mask;
            msb_xb[i] = (data[i] >> shift) & 1;
            tmp_x[i] = (shift_mask - tmp_x[i]) & shift_mask;
        }
    }

    auto millionaire_with_eq = std::make_unique<MillionaireWithEquality>(party, iopack, otpack);
    millionaire_with_eq->compare_with_eq(res_cmp, res_eq, tmp_x.get(), num_cmps, bitlength_data - 1, true, radix_base); // computing greater_than

    for (int i = 0; i < num_cmps; i++)
    {
        res_cmp[i] = res_cmp[i] ^ msb_xb[i] ^ res_eq[i];
    }
    if (greater_than == true && party == sci::ALICE)
    {
        for (int i = 0; i < num_cmps; i++)
        {
            res_cmp[i] = res_cmp[i] ^ 1;
        }
    }
}

void MyFunctions::compare_const_arithmetic(uint8_t *res_cmp, uint64_t *res_cmp_arithmetic, uint64_t *data, uint64_t value, int num_cmps, int bitlength_data, int bitlength_res_arithmetic, bool greater_than, bool equality, int radix_base)
{
    this->compare_const(res_cmp, data, value, num_cmps, bitlength_data, greater_than, equality, radix_base);
    aux->B2A(res_cmp, res_cmp_arithmetic, num_cmps, bitlength_res_arithmetic);
}

void MyFunctions::check_equality_arithmetic(uint8_t *res_eq, uint64_t *res_eq_arithmetic, uint64_t *data, uint64_t value, int num_eqs, int bitlength_data, int bitlength_res_arithmetic, int radix_base)
{
    this->check_equality(res_eq, data, value, num_eqs, bitlength_data, radix_base);
    aux->B2A(res_eq, res_eq_arithmetic, num_eqs, bitlength_res_arithmetic);
}

void MyFunctions::compare_with_eq_const_arithmetic(uint8_t *res_cmp, uint8_t *res_eq, uint64_t *res_cmp_arithmetic, uint64_t *res_eq_arithmetic, uint64_t *data, uint64_t value, int num_cmps, int bitlength_data, int bitlength_res_arithmetic, bool greater_than, int radix_base)
{
    this->compare_with_eq_const(res_cmp, res_eq, data, value, num_cmps, bitlength_data, greater_than, radix_base);

    auto res_cmp_eq = std::make_unique<uint8_t[]>(num_cmps * 2);
    auto res_cmp_eq_arithmetic = std::make_unique<uint64_t[]>(num_cmps * 2);
    memcpy(res_cmp_eq.get(), res_cmp, num_cmps * sizeof(uint8_t));
    memcpy(res_cmp_eq.get() + num_cmps, res_eq, num_cmps * sizeof(uint8_t));
    aux->B2A(res_cmp_eq.get(), res_cmp_eq_arithmetic.get(), num_cmps * 2, bitlength_res_arithmetic);
    memcpy(res_cmp_arithmetic, res_cmp_eq_arithmetic.get(), num_cmps * sizeof(uint64_t));
    memcpy(res_eq_arithmetic, res_cmp_eq_arithmetic.get() + num_cmps, num_cmps * sizeof(uint64_t));
}

void MyFunctions::softmax(int32_t dim, uint64_t *x, uint64_t *y, uint64_t max_x, int32_t bw_x,
                          int32_t bw_y, int32_t s_x, int32_t s_y)
{
    uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
    uint64_t mask_y = (bw_y == 64 ? -1 : ((1ULL << bw_y) - 1));
    uint64_t bw_w = bw_x + 6;
    uint64_t mask_wide = ((bw_w) == 64 ? -1 : ((1ULL << (bw_w)) - 1));
    uint8_t *zero_shares = new uint8_t[dim];

    uint64_t *tmp_1 = new uint64_t[dim];
    uint64_t *tmp_2 = new uint64_t[dim];

    for (int i = 0; i < dim; i++)
    {
        zero_shares[i] = 0;
        x[i] = (x[i] - max_x) & mask_x;
    }
    uint64_t *exp_neg_x = new uint64_t[dim];
    math->lookup_table_exp(dim, x, tmp_2, bw_x, bw_w, s_x, s_y); // exp( neg_x )

    uint64_t sum_exp_x = 0;
    for (int i = 0; i < dim; i++)
    {
        sum_exp_x = ((sum_exp_x + tmp_2[i]) & mask_wide);
    }

    std::cout << "sum_exp_x = " << sum_exp_x << std::endl;

    // sum_exp_x-=8192*dim;
    uint64_t *denominator = new uint64_t[dim];
    for (int i = 0; i < dim; i++)
        denominator[i] = sum_exp_x;

    int scale = s_y + 2;
    int bw_res = s_y + 2;
    uint64_t mask_res = (bw_res == 64 ? -1 : ((1ULL << bw_res) - 1));
    auto res_myscale = std::make_unique<uint64_t[]>(dim);

    math->div(dim, tmp_2, denominator, res_myscale.get(), bw_w, bw_w, bw_res, scale, scale, s_y + 2, false, true);
    std::cout << "bw_nm=" << s_y + 2 << " bw_dn=" << bw_w << " bw_res=" << bw_res << " s_nm=" << scale << " s_dn=" << scale << " s_res=" << s_y << std::endl;
    // for(int i=0;i<dim;i++){
    //     std::cout<<"["<<i<<"]"<<" tmp_2 = "<<tmp_2[i]<<" denominator = " << denominator[i] <<" ans "<< y[i]<<std::endl;
    // }

    math->div(dim, tmp_2, denominator, y, s_y + 2, bw_w, bw_res, s_y, s_y, s_y, false, true);
    std::cout << "bw_nm=" << s_y + 2 << " bw_dn=" << bw_w << " bw_res=" << bw_res << " s_nm=" << s_y << " s_dn=" << s_y << " s_res=" << s_y << std::endl;
    // for(int i=0;i<dim;i++){
    //     std::cout<<"["<<i<<"]"<<" tmp_2 = "<<tmp_2[i]<<" denominator = " << denominator[i] <<" ans "<< y[i]<<std::endl;
    // }

    if (party == sci::ALICE)
    {
        iopackArr[0]->io->send_data(tmp_2, dim * sizeof(uint64_t));
        iopackArr[0]->io->send_data(denominator, dim * sizeof(uint64_t));
        iopackArr[0]->io->send_data(res_myscale.get(), dim * sizeof(uint64_t));
        iopackArr[0]->io->send_data(y, dim * sizeof(uint64_t));
    }
    else if (party == sci::BOB)
    {
        auto recon_nm = std::make_unique<uint64_t[]>(dim);
        auto recon_dn = std::make_unique<uint64_t[]>(dim);
        auto recon_res_myscale = std::make_unique<uint64_t[]>(dim);
        auto recon_res = std::make_unique<uint64_t[]>(dim);
        iopackArr[0]->io->recv_data(recon_nm.get(), dim * sizeof(uint64_t));
        iopackArr[0]->io->recv_data(recon_dn.get(), dim * sizeof(uint64_t));
        iopackArr[0]->io->recv_data(recon_res_myscale.get(), dim * sizeof(uint64_t));
        iopackArr[0]->io->recv_data(recon_res.get(), dim * sizeof(uint64_t));
        for (int i = 0; i < dim; i++)
        {
            recon_nm[i] = (recon_nm[i] + tmp_2[i]) % mask_wide;
            recon_dn[i] = (recon_dn[i] + denominator[i]) % mask_wide;
            recon_res_myscale[i] = (recon_res_myscale[i] + res_myscale[i]) % mask_res;
            recon_res[i] = (recon_res[i] + y[i]) % mask_res;
            std::cout << "[" << i << "]"
                      << " nm = " << recon_nm[i] << " dn = " << recon_dn[i]
                      << " res_myscale = " << recon_res_myscale[i] << " res = " << recon_res[i]
                      << std::endl;
        }
    }

    if (bw_y <= (s_y + 2))
    {
        for (int i = 0; i < dim; i++)
            y[i] = tmp_1[i] & mask_y;
    }
    else
    {
        // xt->z_extend(dim, tmp_1, y, s_y+2, bw_y, zero_shares);
    }

    delete[] tmp_1;
    delete[] tmp_2;
    delete[] exp_neg_x;
}
