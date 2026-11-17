import unittest,json
from pathlib import Path
from variant_cohort.model import analyze,normalize
from variant_cohort.vcf import reference
ROOT=Path(__file__).resolve().parents[1]
def fixture():return json.loads((ROOT/'fixtures/input.json').read_text())
class Variants(unittest.TestCase):
 def test_multiallelic_denominator(self):
  a=analyze(fixture())['alleles'];self.assertEqual(len(a),2)
  self.assertTrue(all(x['allele_number']==5 and x['allele_count']==2 and x['frequency']=='2/5' for x in a));self.assertIsNone(a[0]['dosages']['S3'])
 def test_reference_mismatch(self):
  d=fixture();d['fasta']='>chr1\nAAAAAAAA\n'
  with self.assertRaisesRegex(ValueError,'reference mismatch'):analyze(d)
 def test_genotype_out_of_range(self):
  d=fixture();d['vcf']=d['vcf'].replace('2/2','3/3')
  with self.assertRaisesRegex(ValueError,'out of range'):analyze(d)
 def test_fully_missing(self):
  d=fixture();d['vcf']=d['vcf'].replace('0/1\t2/2\t./1','./.\t.\t.|.')
  for a in analyze(d)['alleles']:self.assertEqual(a['allele_number'],0);self.assertIsNone(a['frequency'])
 def test_haploid_and_phase(self):
  d=fixture();d['vcf']=d['vcf'].replace('0/1\t2/2\t./1','1\t2|0\t0')
  a=analyze(d)['alleles'];self.assertTrue(all(x['allele_number']==4 and x['allele_count']==1 for x in a))
 def test_trim_anchored(self):self.assertEqual(normalize(10,'ACG','ATG'),(11,'C','T'));self.assertEqual(normalize(10,'AA','A'),(10,'AA','A'))
 def test_duplicate_contig(self):
  with self.assertRaises(ValueError):reference('>x\nA\n>x\nC\n')
 def test_duplicate_variant(self):
  d=fixture();d['vcf']+=d['vcf'].splitlines()[-1]+'\n'
  with self.assertRaisesRegex(ValueError,'duplicate variant'):analyze(d)
 def test_empty_records(self):
  d=fixture();d['vcf']='\n'.join(d['vcf'].splitlines()[:2])+'\n';self.assertEqual(analyze(d)['alleles'],[])
 def test_normalized_collision(self):
  d={'fasta':'>x\nACG\n','vcf':'##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nx\t1\t.\tAC\tAT\t.\tPASS\t.\tGT\t0/1\nx\t2\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\n'}
  with self.assertRaisesRegex(ValueError,'collision'):analyze(d)
