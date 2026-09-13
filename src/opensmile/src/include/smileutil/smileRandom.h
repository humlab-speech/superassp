/*  smileRandom.h -- deterministic PRNG for openSMILE components that need noise.
 *
 *  rand()/srand() are forbidden in R packages: the sequence is not portable
 *  between C libraries and the state is not reproducible under R's RNG control.
 *  The components that use randomness (signalGenerator, vectorOperation,
 *  maxIndex) are noise sources, not estimators, so they are given a
 *  self-contained generator instead: MINSTD (Lehmer, a = 16807, m = 2^31-1),
 *  which is exactly what macOS's rand() implements and therefore leaves output
 *  on that platform unchanged.
 */

#ifndef SMILE_RANDOM_H_
#define SMILE_RANDOM_H_

#define SMILE_RANDOM_MAX 2147483647

#ifdef __cplusplus
extern "C" {
#endif

void smile_random_seed(unsigned int seed);
int  smile_random_uniform(void);   /*  in [1, SMILE_RANDOM_MAX - 1], like rand()  */

#ifdef __cplusplus
}
#endif

#endif /* SMILE_RANDOM_H_ */
