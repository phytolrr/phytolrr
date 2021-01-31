import random
import unittest
from tools.motifs import *
from tools.pssm_matrix import *


VALID_AMINO = {'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'V', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y'}


class SimulateMotif(object):
    def __init__(self, offset, score):
        self.offset = offset
        self.score = score


class SimulateMatrix(object):
    def __init__(self):
        self.length = 16
        self.pssm = {}
        for amino in VALID_AMINO:
            self.pssm[amino] = [random.random() for i in range(0, self.length)]


class TestLrrSearch(unittest.TestCase):
    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)

    def test_lrr_search(self):
        test_seq = 'MLLPVKLSPGFPWRWRLSKTLTLFQTLILALLLQPWLSAGVAGYKEDMAALLDFKAYVAPGSSGFLSNWNPNDPDPCSWTGITCDSTKFRVSAIRIAGSDCFKAFSLASHGGCNCNSESVSGNGSSCKIRAQLPGSLGNLTWLRVLSLPFNDLYGEVPREIGSLKLLEELELEGNSLSGNLPMELGLLSSLRVVNLGYNRFTGSIPASLSGCPNLQILNLAGNLLNGTLPDFLGKFSNLRGVFLSFNQLSGSISGDIGNNCEFLEHLHLMGNYFTGGIPSNFGNCSRLRSLLLSSNILDGGIPPDLGRLSALEVLDISRNSLSGHIPSDLGNCKNLSVLILSNEYDPQMFGKEHSLMDSFSGFSSDDKGEFNYFEGGIPDSLANLSMIRIIWAPKATLDGPLPKYWGACKSLQMINLGGNFFNGEFPYSFHECKDMYYFDLSSNKLTGVLTEKLLVPCMGLFNVSGNSLSGDIPKFLETGCPPIPSFVVEKQGQLPSNGLYGQWDYSSVYMSFFACNARSGSSMPYLETDNLPIFHDFSWNNFTGSVPLLPIVAPVRLEMDPLYAFLASGNNISGKLPGYVFETCDILSGMILNLSRNSISGELPEVASNGCISMKQLDVSENRIVGFLPPSFGNLLSLESLDLSRNLLSGQIPMQFGQLKNLRYLSLAGNTLTGGIPSGLAQLPSLQVLELSSNFLTGKIPDGFAGLKNLTSILLDNNKLSGQIPSSFSKMTSLSVFNVSFNNLSGPIPRNVTSVRCDSVLGNPLLLESCHLASQSGPSTEQQGQSGSNTQYAYSPSESVSRKNS'
        # test_seq = 'MKTFSSFFLSVTTLFFFSFFSLSFQASPSQSLYREIHQLISFKDVLPDKNLLPDWSSNKNPCTFDGVTCRDDKVTSIDLSSKPLNVGFSAVSSSLLSLTGLESLFLSNSHINGSVSGFKCSASLTSLDLSRNSLSGPVTTLTSLGSCSGLKFLNVSSNTLDFPGKVSGGLKLNSLEVLDLSANSISGANVVGWVLSDGCGELKHLAISGNKISGDVDVSRCVNLEFLDVSSNNFSTGIPFLGDCSALQHLDISGNKLSGDFSRAISTCTELKLLNISSNQFVGPIPPLPLKSLQYLSLAENKFTGEIPDFLSGACDTLTGLDLSGNHFYGAVPPFFGSCSLLESLALSSNNFSGELPMDTLLKMRGLKVLDLSFNEFSGELPESLTNLSASLLTLDLSSNNFSGPILPNLCQNPKNTLQELYLQNNGFTGKIPPTLSNCSELVSLHLSFNYLSGTIPSSLGSLSKLRDLKLWLNMLEGEIPQELMYVKTLETLILDFNDLTGEIPSGLSNCTNLNWISLSNNRLTGEIPKWIGRLENLAILKLSNNSFSGNIPAELGDCRSLIWLDLNTNLFNGTIPAAMFKQSGKIAANFIAGKRYVYIKNDGMKKECHGAGNLLEFQGIRSEQLNRLSTRNPCNITSRVYGGHTSPTFDNNGSMMFLDMSYNMLSGYIPKEIGSMPYLFILNLGHNDISGSIPDEVGDLRGLNILDLSSNKLDGRIPQAMSALTMLTEIDLSNNNLSGPIPEMGQFETFPPAKFLNNPGLCGYPLPRCDPSNADGYAHHQRSHGRRPA'

        logging.info("Generate matrix...")
        matrix = SimulateMatrix()
        import time
        start = time.time()
        ms = lrr_search(matrix, test_seq)
        end = time.time()
        print(matrix.pssm)
        print(str.format("Find LRR end, sequence length {}, LRR count {}, elapse: {}",
                         len(test_seq), len(ms), end - start))


class TestMotifTools(unittest.TestCase):
    def test_get_highest_score_3_motif(self):
        motifs = [
            SimulateMotif(10, 10.0),
            SimulateMotif(26, 11.0),
            SimulateMotif(42, 15.0)
        ]
        score, motifs = get_highest_score_without_overlay(motifs, 16)
        self.assertLess(score, 36.1)
        self.assertGreater(score, 35.9)
        self.assertEqual(3, len(motifs))

    def test_get_highest_score_1_motif(self):
        motifs = [
            SimulateMotif(10, 10.0)
        ]
        score, motifs = get_highest_score_without_overlay(motifs, 16)
        self.assertLess(score, 10.1)
        self.assertGreater(score, 9.9)
        self.assertEqual(1, len(motifs))

    def test_get_highest_score_overlap_1_time(self):
        motifs = [
            SimulateMotif(10, 10.0),
            SimulateMotif(26, 11.0),
            SimulateMotif(42, 15.0),
            SimulateMotif(57, 14.0)
        ]
        score, motifs = get_highest_score_without_overlay(motifs, 16)
        self.assertLess(score, 36.1)
        self.assertGreater(score, 35.9)
        self.assertEqual(3, len(motifs))
        offsets = set([m.offset for m in motifs])
        self.assertEqual({10, 26, 42}, offsets)

    def test_get_highest_score_overlap_mid_1_time(self):
        motifs = [
            SimulateMotif(10, 10.0),
            SimulateMotif(26, 11.0),
            SimulateMotif(41, 15.0),
            SimulateMotif(57, 14.0)
        ]
        score, motifs = get_highest_score_without_overlay(motifs, 16)
        self.assertLess(score, 39.1)
        self.assertGreater(score, 38.9)
        self.assertEqual(3, len(motifs))
        offsets = set([m.offset for m in motifs])
        self.assertEqual({10, 41, 57}, offsets)

    def test_get_highest_score_overlap_many_times_drop_1(self):
        motifs = [
            SimulateMotif(10, 10.0),
            SimulateMotif(26, 11.0),
            SimulateMotif(41, 15.0),
            SimulateMotif(56, 14.0)
        ]
        score, motifs = get_highest_score_without_overlay(motifs, 16)
        self.assertLess(score, 35.1)
        self.assertGreater(score, 34.9)
        self.assertEqual(3, len(motifs))
        offsets = set([m.offset for m in motifs])
        self.assertEqual({10, 26, 56}, offsets)

    def test_get_highest_score_overlap_many_times_drop_2(self):
        motifs = [
            SimulateMotif(10, 10.0),
            SimulateMotif(26, 11.0),
            SimulateMotif(41, 26.0),
            SimulateMotif(56, 14.0)
        ]
        score, motifs = get_highest_score_without_overlay(motifs, 16)
        self.assertLess(score, 36.1)
        self.assertGreater(score, 35.9)
        self.assertEqual(2, len(motifs))
        offsets = set([m.offset for m in motifs])
        self.assertEqual({10, 41}, offsets)

    def test_get_highest_score_overlap_multiple_2(self):
        motifs = [
            SimulateMotif(10, 10.0),
            SimulateMotif(26, 11.0),
            SimulateMotif(30, 26.0),
            SimulateMotif(35, 12.0),
            SimulateMotif(40, 15.0),
        ]
        score, motifs = get_highest_score_without_overlay(motifs, 16)

        self.assertLess(score, 36.1)
        self.assertGreater(score, 35.9)
        self.assertEqual(2, len(motifs))
        offsets = set([m.offset for m in motifs])
        self.assertEqual({10, 30}, offsets)

    def test_get_highest_score_overlap_multiple_1_4(self):
        motifs = [
            SimulateMotif(10, 10.0),
            SimulateMotif(26, 11.0),
            SimulateMotif(30, 26.0),
            SimulateMotif(35, 12.0),
            SimulateMotif(42, 16.0),
        ]
        score, motifs = get_highest_score_without_overlay(motifs, 16)

        self.assertLess(score, 37.1)
        self.assertGreater(score, 36.9)
        self.assertEqual(3, len(motifs))
        offsets = set([m.offset for m in motifs])
        self.assertEqual({10, 26,42}, offsets)

    def test_get_highest_score_overlap_many_times(self):
        motifs = [
            SimulateMotif(10, 10.0),
            SimulateMotif(26, 11.0),
            SimulateMotif(30, 26.0),
            SimulateMotif(35, 12.0),
            SimulateMotif(42, 16.0),
            SimulateMotif(100, 27.0),
            SimulateMotif(105, 12.0),
            SimulateMotif(116, 17.0),
        ]
        score, motifs = get_highest_score_without_overlay(motifs, 16)

        expect_offsets = {10, 26, 42, 100, 116}
        expect_score = sum([m.score for m in motifs if m.offset in expect_offsets])

        self.assertLess(score, expect_score + 0.01)
        self.assertGreater(score, expect_score - 0.01)
        self.assertEqual(len(expect_offsets), len(motifs))
        offsets = set([m.offset for m in motifs])
        self.assertEqual(expect_offsets, offsets)


if __name__ == '__main__':
    unittest.main()
