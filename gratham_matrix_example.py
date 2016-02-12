#!/usr/bin/env python
import sys
import pprint

def make_grantham_dict(grantham_mat_file):

        """
	Citation:   http://www.ncbi.nlm.nih.gov/pubmed/4843792
	Provenance: http://www.genome.jp/dbget-bin/www_bget?aaindex:GRAR740104
	.	A	R	N	D	C	Q	E	G	H	I	L	K	M	F	P	S	T	W	Y	V
	A	0	112	111	126	195	91	107	60	86	94	96	106	84	113	27	99	58	148	112	64
	R	112	0	86	96	180	43	54	125	29	97	102	26	91	97	103	110	71	101	77	96
	N	111	86	0	23	139	46	42	80	68	149	153	94	142	158	91	46	65	174	143	133
	D	126	96	23	0	154	61	45	94	81	168	172	101	160	177	108	65	85	181	160	152
	C	195	180	139	154	0	154	170	159	174	198	198	202	196	205	169	112	149	215	194	192
	Q	91	43	46	61	154	0	29	87	24	109	113	53	101	116	76	68	42	130	99	96
	E	107	54	42	45	170	29	0	98	40	134	138	56	126	140	93	80	65	152	122	121
	G	60	125	80	94	159	87	98	0	98	135	138	127	127	153	42	56	59	184	147	109
	H	86	29	68	81	174	24	40	98	0	94	99	32	87	100	77	89	47	115	83	84
	I	94	97	149	168	198	109	134	135	94	0	5	102	10	21	95	142	89	61	33	29
	L	96	102	153	172	198	113	138	138	99	5	0	107	15	22	98	145	92	61	36	32
	K	106	26	94	101	202	53	56	127	32	102	107	0	95	102	103	121	78	110	85	97
	M	84	91	142	160	196	101	126	127	87	10	15	95	0	28	87	135	81	67	36	21
	F	113	97	158	177	205	116	140	153	100	21	22	102	28	0	114	155	103	40	22	50
	P	27	103	91	108	169	76	93	42	77	95	98	103	87	114	0	74	38	147	110	68
	S	99	110	46	65	112	68	80	56	89	142	145	121	135	155	74	0	58	177	144	124
	T	58	71	65	85	149	42	65	59	47	89	92	78	81	103	38	58	0	128	92	69
	W	148	101	174	181	215	130	152	184	115	61	61	110	67	40	147	177	128	0	37	88
	Y	112	77	143	160	194	99	122	147	83	33	36	85	36	22	110	144	92	37	0	55
	V	64	96	133	152	192	96	121	109	84	29	32	97	21	50	68	124	69	88	55	0
	"""

	f = open(grantham_mat_file)

	header = f.next().strip().split('\t')
	idx_to_aa = dict(zip(range(0,len(header)), header))


	grantham_dict = {}
	for line in f:
		fields = line.strip().split('\t')
		from_aa = fields[0]

		for idx, score in enumerate(fields):
			if idx == 0:
				continue
			to_aa = idx_to_aa[idx]
			grantham_dict[(from_aa, to_aa)] = score

	return grantham_dict

if __name__ == "__main__":
    grantham_dict = make_grantham_dict(sys.argv[1])

    pp = pprint.PrettyPrinter(indent=2)
    pp.pprint(grantham_dict)
