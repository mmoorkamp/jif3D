#include <cstdlib>
#include "mine.h"
#include "cppmine.h"

using namespace std;


MINE::MINE(double alpha, double c, int est)
{
  char *ret;

  param.alpha = alpha;
  param.c = c;
  param.est = est;
  score = NULL;

  ret = mine_check_parameter(&param);
  if (ret)
    throw ret;
}


MINE::~MINE()
{
  mine_free_score(&score);
}


void MINE::compute_score(double *x, double *y, int n)
{
  prob.x = x;
  prob.y = y;
  prob.n = n;

  mine_free_score(&score);
  score = mine_compute_score(&prob, &param);
  char *ret = "error in mine_compute_score()";
  if (score == NULL)
    throw ret;
}


double MINE::mic()
{
  char *ret = "no score computed";
  if (score == NULL)
    throw ret;

  return mine_mic(score);
}


double MINE::mas()
{
  char *ret = "no score computed";
  if (score == NULL)
    throw ret;

  return mine_mas(score);
}


double MINE::mev()
{
  char *ret = "no score computed";
  if (score == NULL)
    throw ret;

  return mine_mev(score);
}


double MINE::mcn(double eps)
{
  char *ret = "no score computed";
  if (score == NULL)
    throw ret;

  return mine_mcn(score, eps);
}


double MINE::mcn_general()
{
  char *ret = "no score computed";
  if (score == NULL)
    throw ret;

  return mine_mcn_general(score);
}


double MINE::tic(int norm)
{
  char *ret = "no score computed";
  if (score == NULL)
    throw ret;

  return mine_tic(score, norm);
}


double MINE::gmic(double p)
{
  char *ret = "no score computed";
  if (score == NULL)
    throw ret;

  return mine_gmic(score, p);
}
