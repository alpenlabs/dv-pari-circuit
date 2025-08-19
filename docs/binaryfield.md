### Binary Field

Binary Field is a field that has a characteristic of 2, and its size is a power of 2, denoted as GF(2^m). These fields are constructed using polynomials whose coefficients are from GF(2), which contains only the elements 0 and 1. The construction involves selecting an irreducible polynomial M of degree. The field GF(2^m) is then formed by the set of remainders when polynomials are divided by M.

Arithmetic Operations

- Negation: In binary field negation of an element is the same element i.e. x = -x
- Addition and Subtractions: These operations are carry-less and are performed using bitwise XOR. These operations are identical because of the Negation property.
- Squaring: This is a linear operation since (x + y) ^2 = x^2 + y ^2
- Multiplication: This involves a polynomial multiplication followed by a reduction modulo M
- Trace: It is a function that maps an element from the binary field to GF(2).

### Sect233k1 Parameters

Irreducible Polynomial, M = x^233 + x^74 + 1

m = 233

Scalar Field, r = 0x8000000000000000000000000000069d5bb915bcd46efb1ad5f173abdf

a = 0, b = 1

Cofactor h = 4

Generator (0x17232ba853a7e731af129f22ff4149563a419c26bf50a4c9d6eefad6126, 0x1db537dece819b7f70f555a67c427a8cd9bf18aeb9b56e0c11056fae6a3)

### Base Field Multiplication

How do you multiply two polynomials A(x) and B(x) each of degree (m-1) = 232 ? Naive approach would require finding each of the cross terms a_i.b_j (233 x 233 AND ops) then add them suitably to get a 2(m-1) degree polynomial. Karatsuba requires N^1.585 ⇒ 233^1.585 = 5.6k AND gates. 

We employ an FFT-based approach that requires only around 1500 AND gates for a single base field multiplication. An FFT approach evaluates the polynomials at points in a multiplicative subgroup, then multiplies the evaluation pairs, and finally interpolates the product polynomial from the evaluations. With modular reduction, you finally get the product of the polynomial that you’d expect. So let’s dive into each of these.

1. Finding a Multiplicative Subgroup
    
    To do an FFT, we require evaluating polynomials at points in a multiplicative subgroup. Because the product polynomial has degree 2 (m-1) = 464, we require ≥ 464 many evaluation points. Binary fields are known to have multiplicative subgroups of size 2^m -1, so we choose a subgroup generated from GF(2^9) which has 511 many points. This is to say {g^1, g^2,.., g^511} cycles through every non-zero element of the field GF(2^9) without repeating before the end. The generator ‘g’ for this subgroup is chosen to be 2 and this value satisfies the conditions for a generator.
    
2. Polynomial Evaluation
    
    We need to obtain evaluation of polynomial at ≥ 465 points. Evaluation doesn’t cost any AND gates, just XORs. You can do a single polynomial evaluation at a linear time (linear in the degree of polynomial  - 233 coefficients) with Horner’s rule. This algorithm requires a multiplication by constant and an accumulation (addition) done in GF(2^9).
    *Todo: Expand on multiplication by constant*
    
3. Multiplication of Evaluations - Orbit of Points
    
    You could multiply each of the evaluation points in the full domain that will get you ~465 multiplications in GF(2^9). Alternatively you can multiply the evaluations at a selected set of points (around 52) and from there derive the full set of evaluation points with zero AND gates, thus reducing the number of GF(2^9) multiplications from ~465 to ~52.
    
    How this selective multiplication works:
    
    1. Evaluation at square of a domain-point is equal to Square of Evaluation at that point
        
        A(X) = \Sum over i={0, m} (a_i X^i)
        
        A(X^2) = \Sum over i={0, m} (a_i X^2i) 
        
        =  \Sum over i={0, m} (a_i^2 X^2i)  because a_i^2 = a_i in GF{2}
        
        =  (\Sum over i={0, m} (a_i X^i) ) ^2 because of Squaring relation: (x+y)^2 = x^2 + y^2
        
        = A(X)^2
        
        So if you are able to obtain evaluation at a point X, you can obtain evaluation at point X^2 by just squaring the evaluation at X. You repeat this squaring iteratively to obtain evaluations for a set of points {X, X^2, .., X^2k} called an “orbit”.
        
    2. What is an orbit
        
        In GF(2^m), if you iteratively square an element, you form a closed cyclic group, which is called an orbit. Given an element x \in GF(2^m), the set formed by squaring {x, x^2, x^4, x^8, …} is finite because GF(2^m) is finite and eventually x^{2^k} = x, and thus the list closes up into a cycle (orbit) of length dividing m. Two elements x and y lie in the same orbit precisely when y = x^{2^k} for some k.
        
        Sect233k1 has 56 orbits of length 9 and 2 orbits of length 3 and unity (1 orbit of length 1) to give 56*9 + 2*3 + 1 =511 elements of the subgroup. As such, you can have 465 evaluations by starting with evaluation over 52 representatives of each of the orbits of size 9, then square each of these representatives to obtain evaluation over the entire orbit.
        
4. Interpolation
    
    Once you have all evaluations, you interpolate back to obtain a polynomial with degree 2(m-1). This step also requires only XOR gates. A naive version requires O(N^2) where N = 465, to run interpolation. You can lower this cost through alternative approaches.
    Right now we use Karatsuba over Mixed Radix FFT. However the solution is not elegant. We can do much better with Binary Additive IFFT. [Git Issue](https://github.com/alpenlabs/dv-pari-circuit/issues/5)
    
5. Trace and Modular Reduction
    
    Each of the coefficients are in GF(2^9) after the above step. You calculate trace of these elements to obtain values in GF(2) and finally you run modular reduction to obtain values within the range. Both of these require zero AND gates.
    
    *Todo: expand on Trace and Modular reduction*
    

### Scalar Field Multiplication

To multiply scalar field elements of sect233k1, we use Karatsuba to multiply two 232-bit numbers followed by modular reduction that makes efficient use of the structure of field. The modulus is a pseudo-mersenne prime, meaning it’s very close to a power of 2. 

We define C = r - 2^231, where r is modulus and C is a 115-bit constant. 

From this relation we obtain, 2^231 = -C (mod r). This allows us to replace any instance of the large number 2^231 with much smaller 115-bit number C. We use this substitution iteratively until the product of two 232-bit numbers obtained from Karatsuba falls within range.

### Point Addition

1. Tau-Adic Representation
    
    τ-adic representation is a method of expressing an integer as a sum of powers of the complex number τ i.e. an integer n is represented as \Sum over j=[0,l) u_j.τ^j, where the coefficients depend upon scheme used i.e plain {0,1} or naf representations {-1,0,1}.
    
    τ = (\mew + sqrt(-7))/2 where mew=-1 for sect23kk1 and value of l depends upon the algorithm used — with reduced tau-adic representation we get l=235, while we get l = ~460 without any work. The value of ‘l’ is important because it determines how many point additions you will need to do during double-and-add point scalar multiplication algorithm. 
    
    Conversion from a 233-bit big integer to its tau-adic representation is done by using the property “*Lemma28:* *The element c*0 +*c*1 τ *of* Z[τ] *is divisible by* τ *if and only if c*0 *is even*”.
    
    (c0 + c1 τ)/τ = (-c0 + 2 c1)/2 - (c0/2)τ
    
    Generating a PRTNAF of an integer is a bit more complicated but offers significant savings (l = ~235 or 40 digit 5-window NAF representation) 
    
2. Windowed Tau-Adic
    
    With windowed-double-and-add algorithm, you group tau-digits into buckets and use that value as an index for lookup. Windowed version reduces AND gates by w factor at the expense of 2^w precompute table
    
3. Precompute and Lookup Cost
    
    With windowed version you have a precompute table of 2^w rows, which requires around 2^w additions to generate the points in the table — {[0]P, [1]P,[2]P,..[2^w-1]P}. You could get a smaller cost by using point doublings where possible, but the cost won’t differ significantly.
    
    As for lookup, the number of AND gates is (w^2 -1)x233x4. The factor 233x4 comes from how many bits we’re retrieving from the table. Since we represent the co-ordinates in (X,S,Z,T) form — we have 4 233-bit values. Had we used a different representation to reduce the cost say using affine co-ordinate (x,y) or x co-ordinate only, then we’d have to use field inversion somewhere which causes a net increase in AND gates.
    
4. Point frob and Add for Point Scalar Multiplication
    
    Because of the tau-adic representation, we can use frobenius-and-add. Frobenius itself doesn’t cost any AND gate.
    
    ```rust
    pub fn point_scalar_mul(k: &BigUint, point_p: &CurvePoint, w: usize) -> CurvePoint {
        let mut digits = get_tau_adic_digits(&kuint);
        let w_digits = group_tau_adic_digits_into_buckets(&digits, w);
        let table = precompute_lookup_table(point_p, w);
        let mut acc = InnerPoint::identity();
        for &d in &w_digits {
            for _ in 0..w {
                acc = point_frob(&acc);
            }
            let d_mul_p = lookup_table(d)
            acc = point_add(&acc, &d_mul_p);
        }
        acc
    }
    ```
    

### Circuit Builder

1. How circuits are being built
    
    We assign a unique wire-index to each unique wire; logic gates take these labels as reference to the wire. A circuit is just a collection of gates that have reference to some wire index. Gates are defined by mcircuit::Operation structure. We wrap these data on a CktBuilder structure which exposes basic binary operations like xor, and, or. All higher functions are then composed from this primitive functions.
    
    Generating a circuit is akin to compiling a program and as such can’t make assumption about the concrete value that will later be carried by the wire during runtime. 
    
2. Using templates to speed up compilation
    
    Generating a circuit is a sequential process because succeeding gates depend upon wire indices output by preceding gates. This makes compilation process time consuming. A way to speed up this process is through ‘Template’s. Templates are known wire configurations for some time-consuming function. Once a template is generated, we can obtain a unique circuit for the same function by just assigning unique wire indices that still adhere to the wire and gate configuration of the original function — essentially  a clone operation. Clone instead of redo-computation saves 4 times the wait-time because you can use parallel-processing during clone.
    
3. Evaluation and Export to Bristol Fashion
    
    The entire circuit is built from mcircuit::Operation, which supports export to Bristol Fashion and circuit evaluation.


### ~Benchmarks

1. Gate counts for DV Snark Verifier
    - GF9 mul - 33 AND gates, 158 XOR gates
    - GF233 mul - 2124 AND gates, 1,823,131 XOR gates
    - Point Add - 14,868 AND gates, 12,764,014 XOR gates
    - Tau Adic Representation - 438,510 AND gates, 766,570 XOR gates
    - Precompute Table - 2^5 x Point Add
    - Table Lookup - 28892 AND gates, 57,784 XOR gates
    - Point Scalar Mul - 5,012,858 AND gates, 1,602,583,552 billion XOR gates (approx)
    - Field Scalar Mul - 150k AND gates
    - DVSnark Verifier - 5* 2 + 0.15 * 7 = 11M AND gates
2. Compilation Time and Resource Usage
    - Compilation takes around 13 mins
3. Improvements
    - Use PRTNAF to reduce gate count. Preliminary work shows it will reduce the gate count to ~7M from ~10M.
    - Use Binary Additive IFFT for gf233 mul interpolation for a much more elegant and xor-efficient impelemntation. [Git Issue](https://github.com/alpenlabs/dv-pari-circuit/issues/5)