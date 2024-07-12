In this repository, we give the code we use for performance evaluation.


- "SCI" contains a fork of the SCI library. "SCI/PSACEvaluation" contains the codes for the performance evaluation of our 2PC-based protocols. The SCI library involves a SEAL 3.3.2 library, but we do not use it.

- "lattigo_codes" contains the codes to test the performance of the key separation design and the U-test-HE protocol, which are implemented using the lattigo library. 

- "SEAL_codes" contains the codes used to test the performance of the U-test-HE protocol, which is implemented using the Microsoft SEAL 4.1.1 library. 

Since lattigo only supports plaintext modulo to be a prime number, we are unable to set $t$ to a power of 2. But SEAL is certain  to support plaintext modulo a power of 2 (verified in SEAL_codes/U-test-HE). We have plans to implement the MBFV scheme in SEAL to enable more efficient conversion between the BFV scheme and 2PC-based protocols on $Z_{2^l}$. We think it is possible. https://github.com/tuneinsight/lattigo/issues/177#issuecomment-1209884514