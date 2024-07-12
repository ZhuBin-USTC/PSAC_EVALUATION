package main

import (
	"fmt"
	"math/rand"
	"time"

	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

type MyEncoder struct {
	params    bfv.Parameters
	encoder   bfv.Encoder
	evaluator bfv.Evaluator
	print     bool
}

// func (myencoder MyEncoder) bfv_encode_coeffs(coeffs []uint64, pt *bfv.Plaintext) {
// 	ecd := myencoder.encoder

// 	ringT := myencoder.params.RingT()

// 	buffT := ringT.NewPoly()

// 	copy(buffT.Buff, coeffs)

// 	ringT.MulScalar(buffT, pt.scale, buffT)

// 	pt = bfv.NewPlaintext(myencoder.params)
// 	myencoder.encoder.ScaleUp(pt_rt, pt)
// 	return pt
// }

func (myencoder MyEncoder) bfv_encode_coeffs(coeffs []uint64) (pt *bfv.Plaintext) {
	pt_rt := bfv.NewPlaintextRingT(myencoder.params)
	copy(pt_rt.Value.Buff, coeffs)
	pt = bfv.NewPlaintext(myencoder.params)
	myencoder.encoder.ScaleUp(pt_rt, pt)
	return pt
}

func (myencoder MyEncoder) bfv_decode_coeffs(pt *bfv.Plaintext) (coeffs []uint64) {
	pt_rt := bfv.NewPlaintextRingT(myencoder.params)
	myencoder.encoder.ScaleDown(pt, pt_rt)
	coeffs = make([]uint64, myencoder.params.N())
	copy(coeffs, pt_rt.Value.Buff)
	if myencoder.print {
		fmt.Println("coeffs[:100]", coeffs[:100])
	}
	return coeffs
}

func generateRandomUint64Array(K uint64, M int) []uint64 {
	arr := make([]uint64, M)
	for i := 0; i < M; i++ {
		arr[i] = rand.Uint64() % K
	}
	return arr
}

func reverseArray(arr []uint64) []uint64 {
	for i := 0; i < len(arr)/2; i++ {
		j := len(arr) - i - 1
		arr[i], arr[j] = arr[j], arr[i]
	}
	return arr
}

func generateUTestHistXY(K uint64, M int) ([]uint64, []uint64) {
	// test_ctsize1()
	rand.Seed(time.Now().UnixNano())
	HistX := generateRandomUint64Array(K, M)
	HistY := generateRandomUint64Array(K, M)
	return HistX, HistY
}

func (myencoder MyEncoder) generateUTestEncodedHistXY(M int, HistX []uint64, HistY []uint64) (mX *bfv.Plaintext, mY *bfv.Plaintext) {
	InputY := reverseArray(HistY)
	mX = myencoder.bfv_encode_coeffs(HistX)
	mY = myencoder.bfv_encode_coeffs(InputY)
	return mX, mY
}

func (myencoder MyEncoder) generateUTestEncoded_pt_mul_XY(M int) (ptX *bfv.Plaintext, ptY *bfv.Plaintext) {
	coeffX := make([]uint64, M)
	coeffY := make([]uint64, M)
	for i := 0; i < M-1; i++ {
		coeffX[i] = 2
	}
	coeffX[M-1] = 1
	for i := 1; i < M; i++ {
		coeffY[i] = 2
	}
	coeffY[0] = 1
	ptX = myencoder.bfv_encode_coeffs(coeffX)
	ptY = myencoder.bfv_encode_coeffs(coeffY)
	return ptX, ptY
}

func (myencoder MyEncoder) mask_ct(ct *bfv.Ciphertext, index int, M int) *bfv.Ciphertext {
	T := myencoder.params.T()
	N := myencoder.params.N()
	mask_arr := make([]uint64, N)
	for i := 0; i < N; i++ {
		mask_arr[i] = rand.Uint64() % T
	}
	mask_arr[index] = 0
	mask_pt := myencoder.bfv_encode_coeffs(mask_arr)
	masked_ct := myencoder.evaluator.AddNew(ct, mask_pt)
	return masked_ct
}

func test_UTest(P uint64, M int, parametersLiteral bfv.ParametersLiteral) {

	paramDef := parametersLiteral
	// paramDef.T = 4096 * 256
	paramDef.T = 0x3ee0001

	params, err := bfv.NewParametersFromLiteral(paramDef)
	if err != nil {
		panic(err)
	}
	kgen := bfv.NewKeyGenerator(params)

	riderSk, riderPk := kgen.GenKeyPair()
	encoder := bfv.NewEncoder(params)
	decryptor := bfv.NewDecryptor(params, riderSk)
	encryptorRiderPk := bfv.NewEncryptor(params, riderPk)
	encryptorRiderSk := bfv.NewEncryptor(params, riderSk)
	reli := bfv.NewRelinearizationKey(params, 5)
	evaluator := bfv.NewEvaluator(params, rlwe.EvaluationKey{Rlk: reli})
	myencoder := MyEncoder{params, encoder, evaluator, true}

	_, _, _, _ = decryptor, encryptorRiderPk, encryptorRiderSk, evaluator
	HistX, HistY := generateUTestHistXY(P, M)

	mX, mY := myencoder.generateUTestEncodedHistXY(M, HistX, HistY)
	ptX, ptY := myencoder.generateUTestEncoded_pt_mul_XY(M)
	fmt.Println("mX", mX)
	_ = myencoder.bfv_decode_coeffs(mX)
	fmt.Println("mY", mY)
	_ = myencoder.bfv_decode_coeffs(mY)
	fmt.Println("ptX", ptX)
	_ = myencoder.bfv_decode_coeffs(ptX)
	fmt.Println("ptY", ptY)
	_ = myencoder.bfv_decode_coeffs(ptY)

	ctmX := encryptorRiderPk.EncryptNew(mX)
	ctmY := encryptorRiderPk.EncryptNew(mY)
	ctmX_mul_ptX := evaluator.MulNew(ctmX, ptX)
	ctmY_mul_ptY := evaluator.MulNew(ctmY, ptY)

	decX := decryptor.DecryptNew(ctmX)
	fmt.Println("decX")
	_ = myencoder.bfv_decode_coeffs(decX)

	ptmX_mul_ptX := decryptor.DecryptNew(ctmX_mul_ptX)
	ptmY_mul_ptY := decryptor.DecryptNew(ctmY_mul_ptY)
	fmt.Println("ptmX_mul_ptX")
	_ = myencoder.bfv_decode_coeffs(ptmX_mul_ptX)
	fmt.Println("ptmY_mul_ptY")
	_ = myencoder.bfv_decode_coeffs(ptmY_mul_ptY)

	ctuX := evaluator.MulNew(ctmX_mul_ptX, mY)
	ctuY := evaluator.MulNew(ctmY_mul_ptY, mX)
	mask_ctuX := myencoder.mask_ct(ctuX, 2*M-2, M)
	mask_ctuY := myencoder.mask_ct(ctuY, M-1, M)
	ptuX := decryptor.DecryptNew(mask_ctuX)
	ptuY := decryptor.DecryptNew(mask_ctuY)
	uX := myencoder.bfv_decode_coeffs(ptuX)
	fmt.Println("uX[2M-2]", uX[2*M-2])
	uY := myencoder.bfv_decode_coeffs(ptuY)
	fmt.Println("uY[M-1]", uY[M-1])
}

func performance_UTest(P uint64, M int, parametersLiteral bfv.ParametersLiteral) time.Duration {
	paramDef := parametersLiteral
	paramDef.T = 0x3ee0001

	params, err := bfv.NewParametersFromLiteral(paramDef)
	if err != nil {
		panic(err)
	}
	kgen := bfv.NewKeyGenerator(params)

	riderSk, riderPk := kgen.GenKeyPair()
	encoder := bfv.NewEncoder(params)
	decryptor := bfv.NewDecryptor(params, riderSk)
	encryptorRiderPk := bfv.NewEncryptor(params, riderPk)
	encryptorRiderSk := bfv.NewEncryptor(params, riderSk)
	reli := bfv.NewRelinearizationKey(params, 5)
	evaluator := bfv.NewEvaluator(params, rlwe.EvaluationKey{Rlk: reli})
	myencoder := MyEncoder{params, encoder, evaluator, false}

	_, _, _, _ = decryptor, encryptorRiderPk, encryptorRiderSk, evaluator
	HistX, HistY := generateUTestHistXY(P, M)
	mX, mY := myencoder.generateUTestEncodedHistXY(M, HistX, HistY)
	ctmX := encryptorRiderPk.EncryptNew(mX)
	ctmY := encryptorRiderPk.EncryptNew(mY)

	start := time.Now()
	ptX, ptY := myencoder.generateUTestEncoded_pt_mul_XY(M)
	ctmX_mul_ptX := evaluator.MulNew(ctmX, ptX)
	ctmY_mul_ptY := evaluator.MulNew(ctmY, ptY)
	ctuX := evaluator.MulNew(ctmX_mul_ptX, mY)
	ctuY := evaluator.MulNew(ctmY_mul_ptY, mX)
	_ = myencoder.mask_ct(ctuX, 2*M-2, M)
	_ = myencoder.mask_ct(ctuY, M-1, M)
	elapsed := time.Since(start)
	// fmt.Println(elapsed)
	return elapsed
}

func main() {
	// test_UTest(10, 5, bfv.PN13QP218)
	// performance_UTest(10, 5, bfv.PN13QP218)
	time_acc := time.Duration(0)
	for i := 0; i < 10; i++ {
		// time_acc += performance_UTest(30, 100, bfv.PN11QP54)
		// time_acc += performance_UTest(30, 100, bfv.PN13QP218)
		time_acc += performance_UTest(10000, 4000, bfv.PN15QP880)
	}
	fmt.Println(time_acc / 10)
}
