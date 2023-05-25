all : gas gaussian

OUTPUTS_AVERAGED := \
slabcc_CXCHG.dat  slabcc_NXPOT.dat \
slabcc_CXPOT.dat  slabcc_DXCHG.dat  slabcc_MXCHG.dat  slabcc_NYCHG.dat \
slabcc_CYCHG.dat  slabcc_DXPOT.dat  slabcc_MXPOT.dat  slabcc_NYPOT.dat \
slabcc_CYPOT.dat  slabcc_DYCHG.dat  slabcc_MYCHG.dat  slabcc_NZCHG.dat \
slabcc_CZCHG.dat  slabcc_DYPOT.dat  slabcc_MYPOT.dat  slabcc_NZPOT.dat \
slabcc_CZPOT.dat  slabcc_DZCHG.dat  slabcc_MZCHG.dat  slabcc.out       \
slabcc_DIEL.dat   slabcc_NXCHG.dat  slabcc_DZPOT.dat  slabcc_MZPOT.dat

OUTPUTS_RAW := \
slabcc_M.LOCPOT  slabcc_M.CHGCAR slabcc_D.LOCPOT slabcc_D.CHGCAR

% : slabcc.%.in
	../bin/slabcc --input $<
	set -o errexit && for i in $(OUTPUTS_AVERAGED); do numdiff -s' \t\n\\\";' -r 1e-3 -a 1e-5 $$i output/$@/$$i; done
	set -o errexit && for i in $(OUTPUTS_RAW); do numdiff -s' \t\n\\\";' -r 1e-3 -a 1e-4 $$i output/$@/$$i; done

clean:
	rm -f $(OUTPUTS_AVERAGED) $(OUTPUTS_RAW) slabcc.out.* slabcc.log