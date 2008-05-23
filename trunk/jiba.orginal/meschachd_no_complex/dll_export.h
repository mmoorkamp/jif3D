extern	char	*malloc(), *calloc(), *realloc();
extern	void	*malloc(size_t),
extern void m_version();
extern	VEC *v_get(), *v_resize();
extern	MAT *m_get(), *m_resize();
extern	PERM *px_get(), *px_resize();
extern	IVEC *iv_get(), *iv_resize();
extern	int m_free(),v_free();
extern  int px_free();
extern  int iv_free();
extern  BAND *bd_get(), *bd_resize();
extern  int bd_free();
extern	VEC *v_get(int), *v_resize(VEC *,int);
extern	MAT *m_get(int,int), *m_resize(MAT *,int,int);
extern	PERM *px_get(int), *px_resize(PERM *,int);
extern	IVEC *iv_get(int), *iv_resize(IVEC *,int);
extern  BAND *bd_get(int,int,int), *bd_resize(BAND *,int,int,int);
extern  int iv_free(IVEC *);
extern	m_free(MAT *),v_free(VEC *),px_free(PERM *);
extern   int bd_free(BAND *);
extern	void v_foutput(),m_foutput(),px_foutput();
extern  void iv_foutput();
extern	VEC *v_finput();
extern	MAT *m_finput();
extern	PERM *px_finput();
extern	IVEC *iv_finput();
extern	int fy_or_n(), fin_int(), yn_dflt(), skipjunk();
extern	double fin_double();
extern	MAT	*_m_copy(), *m_move(), *vm_move();
extern	VEC	*_v_copy(), *v_move(), *mv_move();
extern	PERM	*px_copy();
extern	IVEC	*iv_copy(), *iv_move();
extern  BAND    *bd_copy();
extern	MAT	*_m_copy(MAT *in,MAT *out,u_int i0,u_int j0),
extern	VEC	*_v_copy(VEC *in,VEC *out,u_int i0),
extern	PERM	*px_copy(PERM *in,PERM *out);
extern	IVEC	*iv_copy(IVEC *in,IVEC *out),
extern  BAND    *bd_copy(BAND *in,BAND *out);
extern	VEC     *v_zero(), *v_rand(), *v_ones();
extern	MAT     *m_zero(), *m_ident(), *m_rand(), *m_ones();
extern	PERM    *px_ident();
extern  IVEC    *iv_zero();
extern	VEC     *v_zero(VEC *), *v_rand(VEC *), *v_ones(VEC *);
extern	MAT     *m_zero(MAT *), *m_ident(MAT *), *m_rand(MAT *),
extern	PERM    *px_ident(PERM *);
extern  IVEC    *iv_zero(IVEC *);
extern	VEC *sv_mlt(), *mv_mlt(), *vm_mlt(), *v_add(), *v_sub(),
extern	double	v_min(), v_max(), v_sum();
extern	VEC	*v_star(), *v_slash(), *v_sort();
extern	double _in_prod(), __ip__();
extern	void	__mltadd__(), __add__(), __sub__(), 
extern	VEC	*sv_mlt(double,VEC *,VEC *),	/* out <- s.x */
extern	double	v_min(VEC *, int *), 
extern	VEC	*v_star(VEC *, VEC *, VEC *),
extern	double	_in_prod(VEC *x,VEC *y,u_int i0), 
extern	void	__mltadd__(Real *,Real *,double,int),
extern	double	_v_norm1(), _v_norm2(), _v_norm_inf(),
extern	double	_v_norm1(VEC *x,VEC *scale),   
extern double m_norm1(MAT *A), m_norm_inf(MAT *A), m_norm_frob(MAT *A);
extern	MAT *sm_mlt(), *m_mlt(), *mmtr_mlt(), *mtrm_mlt(), *m_add(), *m_sub(),
extern   BAND *bd_transp();
extern	MAT *px_rows(), *px_cols(), *swap_rows(), *swap_cols(),
extern	VEC *get_row(), *get_col(), *sub_vec(),
extern	MAT	*sm_mlt(double s,MAT *A,MAT *out), 	/* out <- s.A */
extern  BAND    *bd_transp(BAND *in, BAND *out);   /* out <- A^T */
extern	MAT	*px_rows(PERM *px,MAT *A,MAT *out),	/* out <- P.A */
extern	VEC	*get_row(MAT *,u_int,VEC *),
extern	PERM *px_mlt(), *px_inv(), *px_transp();
extern	int  px_sign();
extern	PERM	*px_mlt(PERM *px1,PERM *px2,PERM *out),	/* out <- px1.px2 */
extern	int	px_sign(PERM *);
extern	IVEC	*iv_add(), *iv_sub(), *iv_sort();
extern	IVEC	*iv_add(IVEC *ix,IVEC *iy,IVEC *out),  /* out <- ix + iy */
extern	double	square(), cube(), mrand();
extern	void	smrand(), mrandlist();
extern  void    m_dump(), px_dump(), v_dump(), iv_dump();
extern MAT *band2mat();
extern BAND *mat2band();
