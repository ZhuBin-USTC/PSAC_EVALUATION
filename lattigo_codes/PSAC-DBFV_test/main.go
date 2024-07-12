package mycode

import (
	"encoding/json"
	"flag"
	"fmt"
	"runtime"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/dbfv"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters). Overrides -short and requires -timeout=0.")
var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")

func testString(opname string, parties int, params bfv.Parameters) string {
	return fmt.Sprintf("%s/LogN=%d/logQ=%d/parties=%d", opname, params.LogN(), params.LogQP(), parties)
}

type testContext struct {
	params bfv.Parameters

	NParties int

	// Polynomial degree
	n int

	// Polynomial contexts
	ringT *ring.Ring
	ringQ *ring.Ring
	ringP *ring.Ring

	encoder bfv.Encoder

	sk0Shards []*rlwe.SecretKey
	sk0       *rlwe.SecretKey

	sk1       *rlwe.SecretKey
	sk1Shards []*rlwe.SecretKey

	pk0 *rlwe.PublicKey
	pk1 *rlwe.PublicKey

	encryptorPk0 bfv.Encryptor
	decryptorSk0 bfv.Decryptor
	decryptorSk1 bfv.Decryptor
	evaluator    bfv.Evaluator

	crs            drlwe.CRS
	uniformSampler *ring.UniformSampler

	sk_RC       *rlwe.SecretKey
	pk_RC       *rlwe.PublicKey
	decryptorRC bfv.Decryptor
}

func gentestContext(params bfv.Parameters, parties int) (tc *testContext, err error) {

	tc = new(testContext)

	tc.params = params

	tc.NParties = parties

	tc.n = params.N()

	tc.ringT = params.RingT()
	tc.ringQ = params.RingQ()
	tc.ringP = params.RingP()

	prng, _ := utils.NewKeyedPRNG([]byte{'t', 'e', 's', 't'})
	tc.crs = prng
	tc.uniformSampler = ring.NewUniformSampler(prng, params.RingQ())

	tc.encoder = bfv.NewEncoder(tc.params)
	tc.evaluator = bfv.NewEvaluator(tc.params, rlwe.EvaluationKey{})

	kgen := bfv.NewKeyGenerator(tc.params)

	// SecretKeys
	tc.sk0Shards = make([]*rlwe.SecretKey, parties)
	// tc.sk1Shards = make([]*rlwe.SecretKey, parties)

	// 生成为0的私钥
	tc.sk0 = bfv.NewSecretKey(tc.params)
	// tc.sk1 = bfv.NewSecretKey(tc.params)

	ringQP, levelQ, levelP := params.RingQP(), params.QCount()-1, params.PCount()-1
	for j := 0; j < parties; j++ {
		tc.sk0Shards[j] = kgen.GenSecretKey()
		// tc.sk1Shards[j] = kgen.GenSecretKey()
		// 各方生成自己的本地的私钥，然后加起来形成全局的私钥
		ringQP.AddLvl(levelQ, levelP, tc.sk0.Value, tc.sk0Shards[j].Value, tc.sk0.Value)
		// ringQP.AddLvl(levelQ, levelP, tc.sk1.Value, tc.sk1Shards[j].Value, tc.sk1.Value)
	}

	// Publickeys
	tc.pk0 = kgen.GenPublicKey(tc.sk0)
	// tc.pk1 = kgen.GenPublicKey(tc.sk1)

	tc.sk_RC = kgen.GenSecretKey()
	tc.pk_RC = kgen.GenPublicKey(tc.sk_RC)

	tc.encryptorPk0 = bfv.NewEncryptor(tc.params, tc.pk0)
	tc.decryptorSk0 = bfv.NewDecryptor(tc.params, tc.sk0)
	// tc.decryptorSk1 = bfv.NewDecryptor(tc.params, tc.sk1)

	tc.decryptorRC = bfv.NewDecryptor(tc.params, tc.sk_RC)
	return
}

func TestDBFVRC(t *testing.T) {

	var err error

	defaultParams := bfv.DefaultParams[:] // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15
	if testing.Short() {
		defaultParams = bfv.DefaultParams[:2] // the short test suite runs for ring degree N=2^12, 2^13
	}
	if *flagLongTest {
		defaultParams = append(defaultParams, bfv.DefaultPostQuantumParams...) // the long test suite runs for all default parameters
	}
	if *flagParamString != "" {
		var jsonParams bfv.ParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			t.Fatal(err)
		}
		defaultParams = []bfv.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, p := range defaultParams {

		var params bfv.Parameters
		if params, err = bfv.NewParametersFromLiteral(p); err != nil {
			t.Fatal(err)
		}

		var tc *testContext
		N := 2
		if tc, err = gentestContext(params, N); err != nil {
			t.Fatal(err)
		}
		for _, testSet := range []func(tc *testContext, t *testing.T){

			// testKeyswitching,
			testPublicKeySwitchingRC,
			testEncToShares,
			// testRefresh,
			// testRefreshAndTransform,
			// testRefreshAndTransformSwitchParams,
			// testMarshalling,
		} {
			testSet(tc, t)
			runtime.GC()
		}

	}
}

func testPublicKeySwitchingRC(tc *testContext, t *testing.T) {

	sk0Shards := tc.sk0Shards
	// pk1 := tc.pk1
	encryptorPk0 := tc.encryptorPk0
	// decryptorSk1 := tc.decryptorSk1

	t.Run(testString("PublicKeySwitchinRC", tc.NParties, tc.params), func(t *testing.T) {

		type Party struct {
			*dbfv.PCKSProtocol
			s     *rlwe.SecretKey
			share *drlwe.PCKSShare
		}

		pcksParties := make([]*Party, tc.NParties)
		// fmt.Println("NParties", tc.NParties)
		for i := 0; i < tc.NParties; i++ {
			p := new(Party)
			p.PCKSProtocol = dbfv.NewPCKSProtocol(tc.params, 6.36)
			p.s = sk0Shards[i]
			p.share = p.AllocateShare()
			pcksParties[i] = p
		}
		P0 := pcksParties[0]

		coeffs, _, ciphertext := newTestVectors(tc, encryptorPk0, t)

		ciphertextSwitched := bfv.NewCiphertext(tc.params, 1)

		for i, p := range pcksParties {
			p.GenShare(p.s, tc.pk_RC, ciphertext.Value[1], p.share)
			if i > 0 {
				P0.AggregateShares(p.share, P0.share, P0.share)
			}
		}

		P0.KeySwitch(ciphertext, P0.share, ciphertextSwitched)
		verifyTestVectors(tc, tc.decryptorRC, coeffs, ciphertextSwitched, t)
	})
}

func testEncToShares(tc *testContext, t *testing.T) {

	coeffs, _, ciphertext := newTestVectors(tc, tc.encryptorPk0, t)

	type Party struct {
		e2s         *dbfv.E2SProtocol
		s2e         *dbfv.S2EProtocol
		sk          *rlwe.SecretKey
		publicShare *drlwe.CKSShare
		secretShare *rlwe.AdditiveShare
	}

	params := tc.params
	P := make([]Party, tc.NParties)

	for i := range P {
		if i == 0 {
			P[i].e2s = dbfv.NewE2SProtocol(params, 3.2)
			P[i].s2e = dbfv.NewS2EProtocol(params, 3.2)
		} else {
			P[i].e2s = P[0].e2s.ShallowCopy()
			P[i].s2e = P[0].s2e.ShallowCopy()
		}

		P[i].sk = tc.sk0Shards[i]
		P[i].publicShare = P[i].e2s.AllocateShare()
		P[i].secretShare = rlwe.NewAdditiveShare(params.Parameters)
	}

	// The E2S protocol is run in all tests, as a setup to the S2E test.
	for i, p := range P {

		p.e2s.GenShare(p.sk, ciphertext.Value[1], p.secretShare, p.publicShare)
		if i > 0 {
			p.e2s.AggregateShares(P[0].publicShare, p.publicShare, P[0].publicShare)
		}
	}

	P[0].e2s.GetShare(P[0].secretShare, P[0].publicShare, ciphertext, P[0].secretShare)

	t.Run(testString("E2SProtocol", tc.NParties, tc.params), func(t *testing.T) {

		rec := rlwe.NewAdditiveShare(params.Parameters)
		for _, p := range P {
			tc.ringT.Add(&rec.Value, &p.secretShare.Value, &rec.Value)
		}

		ptRt := bfv.NewPlaintextRingT(tc.params)
		ptRt.Value.Copy(&rec.Value)

		assert.True(t, utils.EqualSliceUint64(coeffs, tc.encoder.DecodeUintNew(ptRt)))
	})

	crp := P[0].e2s.SampleCRP(params.MaxLevel(), tc.crs)

	t.Run(testString("S2EProtocol", tc.NParties, tc.params), func(t *testing.T) {
		for i, p := range P {
			p.s2e.GenShare(p.sk, crp, p.secretShare, p.publicShare)
			if i > 0 {
				p.s2e.AggregateShares(P[0].publicShare, p.publicShare, P[0].publicShare)
			}
		}

		ctRec := bfv.NewCiphertext(tc.params, 1)
		P[0].s2e.GetEncryption(P[0].publicShare, crp, ctRec)

		verifyTestVectors(tc, tc.decryptorSk0, coeffs, ctRec, t)
	})
}

func newTestVectors(tc *testContext, encryptor bfv.Encryptor, t *testing.T) (coeffs []uint64, plaintext *bfv.Plaintext, ciphertext *bfv.Ciphertext) {

	prng, _ := utils.NewPRNG()
	uniformSampler := ring.NewUniformSampler(prng, tc.ringT)
	coeffsPol := uniformSampler.ReadNew()
	plaintext = bfv.NewPlaintext(tc.params)
	tc.encoder.Encode(coeffsPol.Coeffs[0], plaintext)
	ciphertext = encryptor.EncryptNew(plaintext)
	return coeffsPol.Coeffs[0], plaintext, ciphertext
}

func verifyTestVectors(tc *testContext, decryptor bfv.Decryptor, coeffs []uint64, ciphertext *bfv.Ciphertext, t *testing.T) {
	require.True(t, utils.EqualSliceUint64(coeffs, tc.encoder.DecodeUintNew(decryptor.DecryptNew(ciphertext))))
}
