all : gas gaussian

OUTPUTS := \
slabcc_CXCHG.dat  slabcc_D.LOCPOT   slabcc_M.LOCPOT   slabcc_NXPOT.dat \
slabcc_CXPOT.dat  slabcc_DXCHG.dat  slabcc_MXCHG.dat  slabcc_NYCHG.dat \
slabcc_CYCHG.dat  slabcc_DXPOT.dat  slabcc_MXPOT.dat  slabcc_NYPOT.dat \
slabcc_CYPOT.dat  slabcc_DYCHG.dat  slabcc_MYCHG.dat  slabcc_NZCHG.dat \
slabcc_CZCHG.dat  slabcc_DYPOT.dat  slabcc_MYPOT.dat  slabcc_NZPOT.dat \
slabcc_CZPOT.dat  slabcc_DZCHG.dat  slabcc_MZCHG.dat  slabcc.out       \
slabcc_D.CHGCAR   slabcc_DZPOT.dat  slabcc_MZPOT.dat \
slabcc_DIEL.dat   slabcc_M.CHGCAR   slabcc_NXCHG.dat

% : slabcc.%.in
	../bin/slabcc --input $<
	set -o errexit && for i in $(OUTPUTS); do numdiff $$i output/$@/$$i; done

clean:
	rm -f $(OUTPUTS) slabcc.out.* slabcc.log