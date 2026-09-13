/*  smileRandom.cpp -- see smileRandom.h for the rationale.  */

#include <smileutil/smileRandom.h>

static unsigned int smile_random_state = 1;

void smile_random_seed(unsigned int seed)
{
  /*  MINSTD degenerates for state 0.  */
  smile_random_state = (seed == 0) ? 1u : (seed % (unsigned int)SMILE_RANDOM_MAX);
}

int smile_random_uniform(void)
{
  smile_random_state =
      (unsigned int)(((unsigned long long)smile_random_state * 16807ULL) %
                     (unsigned long long)SMILE_RANDOM_MAX);
  return (int)smile_random_state;
}
