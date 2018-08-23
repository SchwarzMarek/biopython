import unittest

from Bio import SeqIO
from Bio.SeqFeature import ExactPosition, BeforePosition, AfterPosition, SeqFeature, FeatureLocation, WithinPosition, \
    CompoundLocation
from SeqSlicerWithFeatures import slice_sequence_with_features


class SeqSlicerTestsGenuine(unittest.TestCase):
    def setUp(self):
        self.seq = SeqIO.read('SeqSlicer/NC_005816.gb', format='gb')
        # """
        #  source          1..9609
        #                  /organism="Yersinia pestis biovar Microtus str. 91001"
        #                  /mol_type="genomic DNA"
        #                  /strain="91001"
        #                  /db_xref="taxon:229193"
        #                  /plasmid="pPCP1"
        #                  /biovar="Microtus"
        #  repeat_region   1..1954
        #  gene            87..1109
        #                  /locus_tag="YP_pPCP01"
        #                  /db_xref="GeneID:2767718"
        #  CDS             87..1109
        #                  /locus_tag="YP_pPCP01"
        #                  /note="similar to corresponding CDS from previously
        #                  sequenced pPCP plasmid of Yersinia pestis KIM (AF053945)
        #                  and CO92 (AL109969), also many transposase entries for
        #                  insertion sequence IS100 of Yersinia pestis. Contains
        #                  IS21-like element transposase, HTH domain
        #                  (Interpro|IPR007101)"
        #                  /codon_start=1
        #                  /transl_table=11
        #                  /product="putative transposase"
        #                  /protein_id="NP_995567.1"
        #                  /db_xref="GI:45478712"
        #                  /db_xref="GeneID:2767718"
        #                  /translation="MVTFETVMEIKILHKQGMSSRAIARELGISRNTVKRYLQAKSEP
        #                  PKYTPRPAVASLLDEYRDYIRQRIADAHPYKIPATVIAREIRDQGYRGGMTILRAFIR
        #                  SLSVPQEQEPAVRFETEPGRQMQVDWGTMRNGRSPLHVFVAVLGYSRMLYIEFTDNMR
        #                  YDTLETCHRNAFRFFGGVPREVLYDNMKTVVLQRDAYQTGQHRFHPSLWQFGKEMGFS
        #                  PRLCRPFRAQTKGKVERMVQYTRNSFYIPLMTRLRPMGITVDVETANRHGLRWLHDVA
        #                  NQRKHETIQARPCDRWLEEQQSMLALPPEKKEYDVHLDENLVNFDKHPLHHPLSIYDS
        #                  FCRGVA"
        #  misc_feature    87..959
        #                  /locus_tag="YP_pPCP01"
        #                  /note="Transposase and inactivated derivatives [DNA
        #                  replication, recombination, and repair]; Region: COG4584"
        #                  /db_xref="CDD:34222"
        #  misc_feature    <111..209
        #                  /locus_tag="YP_pPCP01"
        #                  /note="Helix-turn-helix domain of Hin and related
        #                  proteins, a family of DNA-binding domains unique to
        #                  bacteria and represented by the Hin protein of Salmonella.
        #                  The basic HTH domain is a simple fold comprised of three
        #                  core helices that form a right-handed...; Region:
        #                  HTH_Hin_like; cl01116"
        #                  /db_xref="CDD:186341"
        #  misc_feature    438..812
        #                  /locus_tag="YP_pPCP01"
        #                  /note="Integrase core domain; Region: rve; cl01316"
        #                  /db_xref="CDD:194099"
        # """

        self.compound = SeqIO.read('SeqSlicer/arab1.gb', format='gb')
        # """
        #  CDS             join(3462..3615,3698..3978,4077..4307,4408..4797,
        #          4876..5028,5141..5332)
        #          /note="containing similarity to NAM-like proteins
        #          gi|3695378"
        #          /codon_start=1
        #          /evidence=not_experimental
        #          /product="T25K16.1"
        #          /protein_id="AAF26460.1"
        #          /db_xref="GI:6715633"
        # """

        self.uniprot_up = SeqIO.read('SeqSlicer/uniprot_up.dat', format='swiss')
        # """
        # FT   DISULFID      7      ?       Interchain (between basic and acidic
        # FT                                chains). {ECO:0000255}.
        # """

    def test_slice(self):
        # check classic slicing
        # ------==============------
        #     |----------------|
        a = slice_sequence_with_features(self.seq, slice(109, 210))
        fea = [f for f in a.features if f.type != 'repeat_region' and 'CDD:186341' in f.qualifiers['db_xref']]
        self.assertEqual(len(fea), 1)
        self.assertIsInstance(fea[0], SeqFeature)
        self.assertIsInstance(fea[0].location, FeatureLocation)
        self.assertIsInstance(fea[0].location.start, BeforePosition)
        self.assertIsInstance(fea[0].location.end, ExactPosition)
        self.assertEqual((int(fea[0].location.start), int(fea[0].location.end)), (1, len(a) - 1))

    def test_slice_removing(self):
        """
        Test if slice removes feature outside location
        # ------=========-------------------
        #                  |----------|
        """
        a = slice_sequence_with_features(self.seq, slice(109, 210))
        fea = [f for f in a.features if f.type != 'repeat_region' and 'CDD:194099' in f.qualifiers['db_xref']]
        self.assertEqual(len(fea), 0)

    def test_slice_before(self):
        """
        Check if slices to feature -> BeforePosition created
        # ------=========-----------
        #          |----------|
        """
        a = slice_sequence_with_features(self.seq, slice(120, 210))
        fea = [f for f in a.features if f.type != 'repeat_region' and 'CDD:186341' in f.qualifiers['db_xref']]
        self.assertEqual(len(fea), 1)
        self.assertIsInstance(fea[0], SeqFeature)
        self.assertIsInstance(fea[0].location, FeatureLocation)
        self.assertIsInstance(fea[0].location.start, BeforePosition)
        self.assertIsInstance(fea[0].location.end, ExactPosition)
        self.assertEqual((int(fea[0].location.start), int(fea[0].location.end)), (0, len(a) - 1))

    def test_slice_after(self):
        """
        Check if slices to feature -> AfterPosition created
        # ----------==========----------
        #       |----------|
        """
        a = slice_sequence_with_features(self.seq, slice(109, 200))
        fea = [f for f in a.features if f.type != 'repeat_region' and 'CDD:186341' in f.qualifiers['db_xref']]
        self.assertEqual(len(fea), 1)
        self.assertIsInstance(fea[0], SeqFeature)
        self.assertIsInstance(fea[0].location, FeatureLocation)
        self.assertIsInstance(fea[0].location.start, BeforePosition)
        self.assertIsInstance(fea[0].location.end, AfterPosition)
        self.assertEqual((int(fea[0].location.start), int(fea[0].location.end)), (1, len(a)))

    def test_slice_inside(self):
        """
        Check if slices to feature -> AfterPosition created
        # ----------==========----------
        #       |----------|
        """
        a = slice_sequence_with_features(self.seq, slice(120, 200))
        fea = [f for f in a.features if f.type != 'repeat_region' and 'CDD:186341' in f.qualifiers['db_xref']]
        self.assertEqual(len(fea), 1)
        self.assertIsInstance(fea[0], SeqFeature)
        self.assertIsInstance(fea[0].location, FeatureLocation)
        self.assertIsInstance(fea[0].location.start, BeforePosition)
        self.assertIsInstance(fea[0].location.end, AfterPosition)
        self.assertEqual((int(fea[0].location.start), int(fea[0].location.end)), (0, len(a)))

    def test_compound_within(self):
        """
        # ----======------======------=======-------======-----
        #   |-----------------------------------------------|
        """
        a = slice_sequence_with_features(self.compound, slice(3400, 5400))
        fea = [f for f in a.features if 'GI:6715633' in f.qualifiers['db_xref']]
        self.assertEqual(len(fea), 1)
        self.assertIsInstance(fea[0], SeqFeature)
        self.assertIsInstance(fea[0].location, CompoundLocation)
        self.assertIsInstance(fea[0].location.start, ExactPosition)
        self.assertIsInstance(fea[0].location.end, ExactPosition)
        self.assertEqual((int(fea[0].location.start), int(fea[0].location.end)), (61, 1932))

    def test_compound_inside1(self):
        """
        # ----======------======------=======-------======-----
        #           |----|
        """
        a = slice_sequence_with_features(self.compound, slice(3620, 3680))
        fea = [f for f in a.features if 'GI:6715633' in f.qualifiers['db_xref']]
        self.assertEqual(len(fea), 0)

    def test_compound_inside2(self):
        """
        # ----======------======------=======-------======-----
        #            |-------------|
        """
        a = slice_sequence_with_features(self.compound, slice(3620, 4050))
        fea = [f for f in a.features if 'GI:6715633' in f.qualifiers['db_xref']]
        self.assertEqual(len(fea), 1)
        self.assertIsInstance(fea[0], SeqFeature)
        self.assertIsInstance(fea[0].location, FeatureLocation)
        self.assertIsInstance(fea[0].location.start, ExactPosition)
        self.assertIsInstance(fea[0].location.end, ExactPosition)
        self.assertEqual((int(fea[0].location.start), int(fea[0].location.end)), (77, len(a) - 72))

    def test_compound_inside3(self):
        """
        # ----======------======------=======-------======-----
        #            |--------------------------|
        """
        a = slice_sequence_with_features(self.compound, slice(3620, 4350))
        fea = [f for f in a.features if 'GI:6715633' in f.qualifiers['db_xref']]
        self.assertEqual(len(fea), 1)
        self.assertIsInstance(fea[0], SeqFeature)
        self.assertIsInstance(fea[0].location, CompoundLocation)
        self.assertIsInstance(fea[0].location.start, ExactPosition)
        self.assertIsInstance(fea[0].location.end, ExactPosition)
        self.assertEqual((int(fea[0].location.start), int(fea[0].location.end)), (77, len(a) - 43))

    def test_compound_slice_inside(self):
        """
        # ----======------======------=======-------======-----
        #                   |-------------|
        """
        a = slice_sequence_with_features(self.compound, slice(3700, 4100))
        fea = [f for f in a.features if 'GI:6715633' in f.qualifiers['db_xref']]
        self.assertEqual(len(fea), 1)
        self.assertIsInstance(fea[0], SeqFeature)
        self.assertIsInstance(fea[0].location, CompoundLocation)
        self.assertIsInstance(fea[0].location.start, BeforePosition)
        self.assertIsInstance(fea[0].location.end, AfterPosition)
        self.assertEqual((int(fea[0].location.start), int(fea[0].location.end)), (0, len(a)))

    def test_uniprot_unknown(self):
        """
        This should raise an error even when doing normal slice well outside the feature location.
        """
        with self.assertRaises(TypeError):
            self.uniprot_up[1:4]


class SeqSlicerTestArtificial(unittest.TestCase):
    def setUp(self):
        self.seq = SeqIO.read('SeqSlicer/gb_artificial.gbk', format='gb')

    def test_slice(self):
        """
             gene            (687.690)..800
                     /gene="within"
        # --------...========---------
        #      |----------------|
        """
        a = slice_sequence_with_features(self.seq, slice(680, 900))
        fea = [f for f in a.features if f.type == 'gene' and 'within1' in f.qualifiers['gene']]
        self.assertEqual(len(fea), 1)
        self.assertIsInstance(fea[0], SeqFeature)
        self.assertIsInstance(fea[0].location, FeatureLocation)
        self.assertIsInstance(fea[0].location.start, WithinPosition)
        self.assertIsInstance(fea[0].location.end, ExactPosition)
        self.assertEqual((int(fea[0].location.start), int(fea[0].location.end)), (6, 120))

    def test_slice_inner(self):
        """
             gene            (687.690)..800
                     /gene="within"
        # --------...========---------
        #             |----|
        """
        a = slice_sequence_with_features(self.seq, slice(700, 750))
        fea = [f for f in a.features if f.type == 'gene' and 'within1' in f.qualifiers['gene']]
        self.assertEqual(len(fea), 1)
        self.assertIsInstance(fea[0], SeqFeature)
        self.assertIsInstance(fea[0].location, FeatureLocation)
        self.assertIsInstance(fea[0].location.start, BeforePosition)
        self.assertIsInstance(fea[0].location.end, AfterPosition)
        self.assertEqual((int(fea[0].location.start), int(fea[0].location.end)), (0, 50))

    def test_slice_inside_left(self):
        """
             gene            (687.690)..800
                     /gene="within"
        # --------...========---------
        #       |------|
        """
        a = slice_sequence_with_features(self.seq, slice(600, 700))
        fea = [f for f in a.features if f.type == 'gene' and 'within1' in f.qualifiers['gene']]
        self.assertEqual(len(fea), 1)
        self.assertIsInstance(fea[0], SeqFeature)
        self.assertIsInstance(fea[0].location, FeatureLocation)
        self.assertIsInstance(fea[0].location.start, WithinPosition)
        self.assertIsInstance(fea[0].location.end, AfterPosition)
        self.assertEqual((int(fea[0].location.start), int(fea[0].location.end)), (86, 100))

    def test_slice_inside_right(self):
        """
             gene            (687.690)..800
                     /gene="within"
        # --------...========---------
        #               |------|
        """
        a = slice_sequence_with_features(self.seq, slice(700, 900))
        fea = [f for f in a.features if f.type == 'gene' and 'within1' in f.qualifiers['gene']]
        self.assertEqual(len(fea), 1)
        self.assertIsInstance(fea[0], SeqFeature)
        self.assertIsInstance(fea[0].location, FeatureLocation)
        self.assertIsInstance(fea[0].location.start, BeforePosition)
        self.assertIsInstance(fea[0].location.end, ExactPosition)
        self.assertEqual((int(fea[0].location.start), int(fea[0].location.end)), (0, 100))

    def test_slice_to_uncertain(self):
        """
             gene            (687.690)..800
                     /gene="within"
        # --------...========---------
        #          |------|
        """
        a = slice_sequence_with_features(self.seq, slice(688, 750))
        fea = [f for f in a.features if f.type == 'gene' and 'within1' in f.qualifiers['gene']]
        self.assertEqual(len(fea), 0)

    def test_slice_to_uncertain_keep(self):
        """
             gene            (687.690)..800
                     /gene="within"
        # --------...========---------
        #          |------|
        """
        a = slice_sequence_with_features(self.seq, slice(688, 750), keep_all_features=True)
        fea = [f for f in a.features if f.type == 'gene' and 'within1' in f.qualifiers['gene']]
        self.assertEqual(len(fea), 1)
        self.assertIsInstance(fea[0], SeqFeature)
        self.assertIsInstance(fea[0].location, FeatureLocation)
        self.assertIsInstance(fea[0].location.start, WithinPosition)
        self.assertIsInstance(fea[0].location.end, AfterPosition)
        self.assertEqual((int(fea[0].location.start), int(fea[0].location.end)), (0, len(a)))

    def test_slice_to_uncertain_r(self):
        """
        gene            900..(950.960)
                        /gene="within2"
        # --------========...---------
        #           |------|
        """
        a = slice_sequence_with_features(self.seq, slice(920, 955))
        fea = [f for f in a.features if f.type == 'gene' and 'within2' in f.qualifiers['gene']]
        self.assertEqual(len(fea), 0)

    def test_slice_to_uncertain_keep_r(self):
        """
        gene            900..(950.960)
                        /gene="within2"
        # --------========...---------
        #           |------|
        """
        a = slice_sequence_with_features(self.seq, slice(920, 955), keep_all_features=True)
        fea = [f for f in a.features if f.type == 'gene' and 'within2' in f.qualifiers['gene']]
        self.assertEqual(len(fea), 1)
        self.assertIsInstance(fea[0], SeqFeature)
        self.assertIsInstance(fea[0].location, FeatureLocation)
        self.assertIsInstance(fea[0].location.start, BeforePosition)
        self.assertIsInstance(fea[0].location.end, WithinPosition)
        self.assertEqual((int(fea[0].location.start), int(fea[0].location.end)), (0, len(a)))

if __name__ == "__main__":
    unittest.main()