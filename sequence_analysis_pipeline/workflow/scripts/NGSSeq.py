import nwalign as nw
# https://pypi.org/project/nwalign/


class NGSSeq:
    """Class encapsulating information sequences in database"""

    def __init__(self):
        self.database = None
        self.naseq = None
        self.seq = None
        self.length = None
        # segextents; #(:,i,j) - segment i (0:10), j=0 is start, j=1 is end,  end==start-1 for missing segment
        self.segextents = None
        # Class of each sequence(index into seqclasses)
        self.seqclass = seqclass
        self.seqnames = seqnames  # Mapping between naseq and names of sequences
        self.tags = tags  # Tags
        self.hits = hits  # Hits
        self.clusters = clusters  # Master data for each cluster

        # TODO proporties (Transient)
        # (:,i) Identity of each of i=0:10 segments (pre,s3a,s1a,loop1,s1b,core,s2a,loop2,s2b,s3b,post), automatically rebuilt as needed
        self.segments = None
        self.otherseqs = None  # Lookup for temporary seqs read from database
        # Distance between cluster roots (struct:  clusters(N),dists(N,N) )
        self.clusterdists = None

        # TODO proporties (Constant)
        # Cell vector becomes a tuple (https://www.mathworks.com/help/matlab/matlab_external/passing-data-to-python.html)
        # Seq classes
        # Each entry consists of name, total lengt hrange, l1 len range or regexp, l2 len range or regexp, planned?
        # Evaluated in order
        self.seqclasses = (('Short', [0, 40], [], [], False),
                           ('NsN30', [], [1, 9], [28, 32], True),
                           ('N0N30', [], [0, 0], [28, 32], False),
                           ('N30N30', [], [28, 32], [28, 32], False),
                           ('N30Ns', [], [28, 32], [1, 9], True),
                           ('N30N0', [], [28, 32], [0, 0], False),
                           ('NsNs', [], [0, 9], [0, 9], False),
                           # Minimum length of any structured library is 84nt
                           # ('3wj',[84,999],'NNN RYRYRYRYRY NNNNN RYRYRYRYRY NNN RYRYRYRYRY NNNNN RYRYRYRYRY NNN',[],True),
                           # ('4wj',[84,999],'NN RYRYRYRY NNNN RYRYRYRY NN RYRYRYRYRY NNNN RYRYRYRYRY NN RYRYRYRY NNNN RYRYRYRY NN',[],True),
                           # ('5wj',[84,999],'NN RYRYRYRY NNNN RYRYRYRY NN RYRYRYRY NNNN RYRYRYRY NN RYRYRYRY NNNN RYRYRYRY NN RYRYRYRY NNNN RYRYRYRY NN',[],True),
                           # ('GR',[84,999],'NNNNNNN CGCGTGGATATGGCACGCA NNNNNNNNNN GGGCACCGTAAATGTCC NNNNNN',[],True),
                           # ('GRP1d3',[84,999],'NNNN CGCGTGGATATGGCACGCA NNNNNNNNNN GGGCACCGTAAATGTCC NNN',[],True),
                           # ('CG',[84,999],'NNNNNNNNN CAAACCATTCGAAAGAGTGGGACG NNNNN CCTCCGGCCTAAACCGAAAGGTAGGTAGCGGGG NNNNNNN',[],True),
                           # ('CG_P1d2',[84,999],'NNNNNNN CAAACCATTCGAAAGAGTGGGACG NNNNN CCTCCGGCCTAAACCGAAAGGTAGGTAGCGGGG NNNN',[],True),
                           # ('HH',[84,999],'NNNN TGGTATCCAATGAAAATGTACTACCA NNNNNNNNNN CCCAAATAGG NNNNNNN',[],True),
                           # ('HH_P1d2',[84,999],'NN TGGTATCCAATGAAAATGTACTACCA NNNNNNNNNN CCCAAATAGG NNNNN',[],True),
                           # ('GLY',[84,999],'NNNN GAACCGTTTAATCGGTCGC NNNNNNNN CAAGCTCTGCGCATATGCAGAGTGAAA NNNNNNN GCAAAA NNN',[],True),
                           # ('TET',[84,999],'AAAA NNNNN CAGATTTCGATCTG NNNN GGTG NNNNNNNNNNN CACCT',[],True),
                           # ('GLN',[84,999],'N?N?NNN TTGGCCCAGTTTATCTGGG NNNNNNNN AGGTCTTTGGCCT NNNN CAA NNN',[],True),
                           # ('FMN_p1',[84,999],'NNNNCAGGGTGAAATTCCCGANNNNTGGTATAGTCCANNNAAGTATTTGCTTNNNTTTGGTGAAATTCCAAAACNNNNDGTAGAGTCHNNNNNNNNN',[],True),
                           # ('FMN_p4',[84,999],'NNNTTTGGTGAAATTCCAAAACNNNNDGTAGAGTCHNNNNNNNNNGAAGAGAAATCTTCNNNNCAGGGTGAAATTCCCGANNNNTGGTATAGTCCANNN',[],True),
                           # ('SAM_p1',[84,999],'NNNNNNGGTGGAGGGACTGGCCCGATGAAACCCNNNNNNNAGAAATNNNNNNAATTCCTGCAGCGGTTTCGCTGAAA',[],True),
                           # ('SAM_p2',[84,999],'NNNNAATTCCTGCAGCGGTTTCGCTGAAACTNNNNNNNNNN          GGTGGAGGGACTGGCCCGATGAAACCCNNNNN',[],True),
                           ('NsN60', [], [1, 9], [58, 62], True),
                           ('N0N60', [], [0, 0], [58, 62], False),
                           ('N60Ns', [], [58, 62], [1, 9], True),
                           ('N60N0', [], [58, 62], [0, 0], False),
                           ('N60N60', [], [58, 62], [58, 62], False),
                           ('NsNm', [], [0, 9], [10, 29], False),
                           ('NsNl', [], [0, 9], [30, 999], False),
                           ('NmNs', [], [10, 29], [0, 9], False),
                           ('NmNm', [], [10, 29], [10, 29], False),
                           ('NmNl', [], [10, 29], [30, 999], False),
                           ('NlNs', [], [30, 999], [0, 9], False),
                           ('NlNm', [], [30, 999], [10, 29], False),
                           ('NlNl', [], [30, 999], [30, 999], False),
                           ('Unclassifiable', [], [0, 100], [0, 100], False),
                           ('Unparseable', [], [], [], False))

        # Ribozyme parts (empty sets are wildcards)
        self.pts = ((),
                    ('GCTGTC', 'ACTGTC', 'TCTGTC', 'GCTTTC', 'GCTTTTC', 'GCTGTA', 'GCTATC', 'GTTGTC', 'GATGTC', 'GCAGTC',
                     'GCCGTC', 'GCTGTT', 'GCTGCC', 'GCTGGC', 'GCTGAC', 'GCTCTC', 'GCGGTC', 'GCTGTG', 'GGTGTC', 'CTGTC', 'TGTC'),
                    ('ACCGGA', 'ACAGGA', 'ACTGGA', 'AACGAA', 'ACGGGA',
                     'ACCGAA', 'ATCGGA', 'ACCAGA', 'ACAGAA'),
                    (),
                    ('TCCGGT', 'TCCAGT', 'TCTGGT', 'TCCTGT', 'TCCGAT', 'TTCTGT', 'TCTTGT',
                     'TCTAGT', 'TCCCGT', 'TACTGT', 'TCATGT', 'TTCGGT', 'TCTGT'),
                    ('CTGATGA', 'CTGACGA', 'CCGATGA', 'CTGAAGA', 'CTGTCTGA', 'CTCATGA', 'TTGAAGA', 'CTGAGGA', 'CTGATTA', 'CTAATGA', 'ATGATGA', 'CTTATGA', 'CTGATGG',
                     'CGGATGA', 'CTGATGT', 'CTGATAA', 'TTGATGA', 'CTGGTGA', 'TTGATAA', 'CTGTTGA', 'CTAACGA', 'CAGATGA', 'CCGACGA', 'GTGATGA', 'CTGCTGA', 'TGATGA', 'TGACGA'),
                    ('GTCC', 'GTCT', 'ATCC', 'TTCC', 'GTTC', 'GCCC'),
                    (),
                    ('GGAC', 'CGAC'),
                    ('GAAACAGC'),
                    ())
        self.scoringmatrix = NGSSeq.editdistance()

        # TODO fix this method, since struct is not a Python thing
        @staticmethod
        def load_obj(obj):
            """Handler for loading from file"""
            if isstruct(obj):
                print("WarningL loading a struct -- run NGSSeq::cleanup()\n")
            else:
                obj.cleanup()
            return obj

        @staticmethod
        def editdistance():
            """Scoring matrix that counts edit distance
            Use with gapopen=1
            """
            # TODO import nuc44 scoring matrix
            sm = nuc44() > 0
            sm[-1, :] = 1
            sm[:, -1] = 1
            sm = sm - 1

            return sm

        @staticmethod
        def seqcompare(seq1, seq2):
            """Align 2 sequences using scoring based on edit differences
            Substitution count as 1, insertion/deletions as 1"""
            # al stands for alignment
            alseq1, alseq2 = nw.global_align(
                seq1, seq2, matrix=self.scoringmatrix, gap_open=1)
            dist, al = nw.score_alignment(
                alseq1, alseq2, matrix=self.scoringmatrix, gap_open=1)
            dist = -dist
            return dist, al

        def

# # Class encapsulating information sequences in database
# classdef NGSSeq  < matlab.mixin.Copyable
#   properties
#     database;   # Which database this represents (e.g. NGS34)
#     naseq;
#     seq;
#     length;
#     segextents; # (:,i,j) - segment i (1:11), j=1 is start, j=2 is end,  end==start-1 for missing segment
#     class; 	# Class of each sequence (index into seqclasses)
#     seqnames;	# Mapping between naseq and names of sequences
#     tags;	# Tags
#     hits;	# Hits
#     clusters;   # Master data for each cluster
#   end

#   properties (Transient)
#     segments;   # (:,i) Identity of each of i=1:11 segments (pre,s3a,s1a,loop1,s1b,core,s2a,loop2,s2b,s3b,post), automatically rebuilt as needed
#     otherseqs;  # Lookup for temporary seqs read from database
#     clusterdists;   # Distance between cluster roots (struct:  clusters(N),dists(N,N) )
#   end

#   properties (Constant)
#     # Seq classes
#     # Each entry consists of name, totallengthrange, l1lenrange or regexp, l2lenrange or regexp, planned?
#     # Evaluated in order
#     seqclasses=(('Short',[0,40],[],[],False),
#                 ('NsN30',[],[1,9],[28,32],True),
#                 ('N0N30',[],[0,0],[28,32],False),
#                 ('N30N30',[],[28,32],[28,32],False),
#                 ('N30Ns',[],[28,32],[1,9],True),
#                 ('N30N0',[],[28,32],[0,0],False),
#                 ('NsNs',[],[0,9],[0,9],False),
#     # Minimum length of any structured library is 84nt
#                 # ('3wj',[84,999],'NNN RYRYRYRYRY NNNNN RYRYRYRYRY NNN RYRYRYRYRY NNNNN RYRYRYRYRY NNN',[],True),
#                 # ('4wj',[84,999],'NN RYRYRYRY NNNN RYRYRYRY NN RYRYRYRYRY NNNN RYRYRYRYRY NN RYRYRYRY NNNN RYRYRYRY NN',[],True),
#                 # ('5wj',[84,999],'NN RYRYRYRY NNNN RYRYRYRY NN RYRYRYRY NNNN RYRYRYRY NN RYRYRYRY NNNN RYRYRYRY NN RYRYRYRY NNNN RYRYRYRY NN',[],True),
#                 # ('GR',[84,999],'NNNNNNN CGCGTGGATATGGCACGCA NNNNNNNNNN GGGCACCGTAAATGTCC NNNNNN',[],True),
#                 # ('GRP1d3',[84,999],'NNNN CGCGTGGATATGGCACGCA NNNNNNNNNN GGGCACCGTAAATGTCC NNN',[],True),
#                 # ('CG',[84,999],'NNNNNNNNN CAAACCATTCGAAAGAGTGGGACG NNNNN CCTCCGGCCTAAACCGAAAGGTAGGTAGCGGGG NNNNNNN',[],True),
#                 # ('CG_P1d2',[84,999],'NNNNNNN CAAACCATTCGAAAGAGTGGGACG NNNNN CCTCCGGCCTAAACCGAAAGGTAGGTAGCGGGG NNNN',[],True),
#                 # ('HH',[84,999],'NNNN TGGTATCCAATGAAAATGTACTACCA NNNNNNNNNN CCCAAATAGG NNNNNNN',[],True),
#                 # ('HH_P1d2',[84,999],'NN TGGTATCCAATGAAAATGTACTACCA NNNNNNNNNN CCCAAATAGG NNNNN',[],True),
#                 # ('GLY',[84,999],'NNNN GAACCGTTTAATCGGTCGC NNNNNNNN CAAGCTCTGCGCATATGCAGAGTGAAA NNNNNNN GCAAAA NNN',[],True),
#                 # ('TET',[84,999],'AAAA NNNNN CAGATTTCGATCTG NNNN GGTG NNNNNNNNNNN CACCT',[],True),
#                 # ('GLN',[84,999],'N?N?NNN TTGGCCCAGTTTATCTGGG NNNNNNNN AGGTCTTTGGCCT NNNN CAA NNN',[],True),
#                 # ('FMN_p1',[84,999],'NNNNCAGGGTGAAATTCCCGANNNNTGGTATAGTCCANNNAAGTATTTGCTTNNNTTTGGTGAAATTCCAAAACNNNNDGTAGAGTCHNNNNNNNNN',[],True),
#                 # ('FMN_p4',[84,999],'NNNTTTGGTGAAATTCCAAAACNNNNDGTAGAGTCHNNNNNNNNNGAAGAGAAATCTTCNNNNCAGGGTGAAATTCCCGANNNNTGGTATAGTCCANNN',[],True),
#                 # ('SAM_p1',[84,999],'NNNNNNGGTGGAGGGACTGGCCCGATGAAACCCNNNNNNNAGAAATNNNNNNAATTCCTGCAGCGGTTTCGCTGAAA',[],True),
#                 # ('SAM_p2',[84,999],'NNNNAATTCCTGCAGCGGTTTCGCTGAAACTNNNNNNNNNN          GGTGGAGGGACTGGCCCGATGAAACCCNNNNN',[],True),
#                 ('NsN60',[],[1,9],[58,62],True),
#                 ('N0N60',[],[0,0],[58,62],False),
#                 ('N60Ns',[],[58,62],[1,9],True),
#                 ('N60N0',[],[58,62],[0,0],False),
#                 ('N60N60',[],[58,62],[58,62],False),
#                 ('NsNm',[],[0,9],[10,29],False),
#                 ('NsNl',[],[0,9],[30,999],False),
#                 ('NmNs',[],[10,29],[0,9],False),
#                 ('NmNm',[],[10,29],[10,29],False),
#                 ('NmNl',[],[10,29],[30,999],False),
#                 ('NlNs',[],[30,999],[0,9],False),
#                 ('NlNm',[],[30,999],[10,29],False),
#                 ('NlNl',[],[30,999],[30,999],False),
#                 ('Unclassifiable',[],[0,100],[0,100],False),
#                 ('Unparseable',[],[],[],False));

#     # Ribozyme parts (empty sets are wildcards)
#     pts=((),
#          ('GCTGTC','ACTGTC','TCTGTC','GCTTTC','GCTTTTC','GCTGTA','GCTATC','GTTGTC','GATGTC','GCAGTC','GCCGTC','GCTGTT','GCTGCC','GCTGGC','GCTGAC','GCTCTC','GCGGTC','GCTGTG','GGTGTC','CTGTC','TGTC'),
#          ('ACCGGA','ACAGGA','ACTGGA','AACGAA','ACGGGA','ACCGAA','ATCGGA','ACCAGA','ACAGAA'),
#          (),
#          ('TCCGGT','TCCAGT','TCTGGT','TCCTGT','TCCGAT','TTCTGT','TCTTGT','TCTAGT','TCCCGT','TACTGT','TCATGT','TTCGGT','TCTGT'),
#          ('CTGATGA','CTGACGA','CCGATGA','CTGAAGA','CTGTCTGA','CTCATGA','TTGAAGA','CTGAGGA','CTGATTA','CTAATGA','ATGATGA','CTTATGA','CTGATGG','CGGATGA','CTGATGT','CTGATAA','TTGATGA','CTGGTGA','TTGATAA','CTGTTGA','CTAACGA','CAGATGA','CCGACGA','GTGATGA','CTGCTGA','TGATGA','TGACGA'),
#          ('GTCC','GTCT','ATCC','TTCC','GTTC','GCCC'),
#          (),
#          ('GGAC','CGAC'),
#          ('GAAACAGC'),
#          ());
#     scoringmatrix=NGSSeq.editdistance();
#   end
