#ifndef RNG_INCLUDED
#define RNG_INCLUDED

#include <cstdint>

/* xoshiro256++ 1.0, by Blackman & Vigna (public domain).
   Satisfies std::UniformRandomBitGenerator, so it is a drop-in replacement for
   std::mt19937 in the standard distributions.  Chosen over mt19937 because the
   simulator draws two or three variates per attempted move, and profiling showed
   mt19937 accounting for over a quarter of all retired instructions. */
struct Rng {
  typedef uint32_t result_type;
  uint64_t s[4];
  uint32_t spare;
  bool haveSpare;

  static constexpr result_type min() { return 0; }
  static constexpr result_type max() { return 0xffffffffu; }

  explicit Rng (uint64_t seed = 1) : spare(0), haveSpare(false) {
    for (int i = 0; i < 4; ++i) {           /* splitmix64 seeding */
      seed += 0x9e3779b97f4a7c15ULL;
      uint64_t z = seed;
      z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
      z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
      s[i] = z ^ (z >> 31);
    }
    for (int i = 0; i < 8; ++i) next64();   /* discard a little to mix the state */
  }

  static inline uint64_t rotl (uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

  inline uint64_t next64() {
    const uint64_t r = rotl (s[0] + s[3], 23) + s[0];
    const uint64_t t = s[1] << 17;
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl (s[3], 45);
    return r;
  }

  inline result_type operator()() {
    if (haveSpare) { haveSpare = false; return spare; }
    const uint64_t v = next64();
    spare = (uint32_t) (v >> 32);
    haveSpare = true;
    return (uint32_t) v;
  }

  /* unbiased-enough index in [0,n) without a division: Lemire multiply-shift */
  inline uint32_t index (uint32_t n) {
    return (uint32_t) (((uint64_t) (*this)() * (uint64_t) n) >> 32);
  }
  /* uniform double in [0,1) */
  inline double real() { return (next64() >> 11) * 0x1.0p-53; }
};

#endif /* RNG_INCLUDED */
