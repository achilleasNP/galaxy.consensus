import vcf as pyvcf
import itertools as it
from variant_ensemble import variant_ensemble
from functools import total_ordering
from heapq import merge



@total_ordering
class simple_variant:
  '''
  Minimum representation of a variant.
  Dictionary order with letter based chromosomes
  being bigger than number based ones. Number based
  chromosomes compared as integers and alpha ones
  compaired lexicographically.
  '''

  def __init__(self, rec, ALT):
    '''
    Take a pyVCF _Record and make a string. Must provide explicit alternate allele ALT.
    '''
    self.ID = rec.CHROM + ':' + str(rec.POS) + ':' + str(rec.REF) + ':' + str(ALT)
    self.contig = rec.CHROM
    self.pos = rec.POS
    self.alt = ALT 
    self.ref = rec.REF 
  def __hash__(self):
    return hash(self.ID)
  def __eq__(self, other):
    if self.ID == other.ID: return True
    else: return False
  def __str__(self):
    return self.ID
  def __gt__(self,other):
      if self.contig.isalpha() and not other.contig.isalpha():
          out = True
      elif self.contig.isalpha() and other.contig.isalpha():
          out = (self.contig,self.pos,self.ref,self.alt) > (other.contig, other.pos,self.ref,self.alt)
      elif not self.contig.isalpha() and not other.contig.isalpha():
          out = (int(self.contig),self.pos,self.ref,self.alt) > (int(other.contig), other.pos,self.ref,self.alt)
      else:
          out = False
      return out


class vcf_ensemble:
  '''
  Represents an arbitrary number of VCF files describing the same data set.
  '''
  @property
  def samples(self):
    '''
    Samples common to all input VCF files.
    '''
    sampleSets = [ set(x.samples) for x in self.vcfReaders ]
    self.samples = reduce( lambda x,y: x.intersection(y), sampleSets )
    return self.samples

  @samples.setter
  def samples(self, samples):
    self.samples = samples

  @samples.getter
  def samples(self):
    sampleSets = [ set(x.samples) for x in self.vcfReaders ]
    self.samples = reduce( lambda x,y: x.intersection(y), sampleSets )
    return self.samples
  
  @property
  def variants(self):
    '''
    List of variant sets found in each caller.
    '''
    variantSets = []
    readers = [ pyvcf.Reader(open(x), prepend_chr = False) for x in self.vcfs]
    variants = merge(*[get_variant_iterator(reader) for reader in readers ])
    return variants

  @property
  def vcfReaders(self):
    '''
    Connections to VCF files.
    '''
    return [ pyvcf.Reader(open(x), prepend_chr = False) for x in self.vcfs ]


  def __init__(self, *args, **kwargs):
    self.vcfs = kwargs.get('vcfList')
    self.ignoreMissing = kwargs.get('ignoreMissing')

  def _concordant_sites(self, threshold):
    '''
    Find variants that agree at some threshold.

    This problem reduces to finding the union of a set of sets. The top level
    set consists of combinations of VCF sets set by the threshold (i.e.
   variants choose threshold).
    '''
    for _,group in it.groupby(self.variants):
        n_vars = 0
        for var in group:
            n_vars += 1
            if n_vars >= threshold:
                yield var
                break

    
             

  def concordant_variants(self, siteThresh, genoThresh):
    '''
    Reduce the ensemble of VCFs into a set of variant sites that agree at some threshold.
    '''
    variants = self._concordant_sites(siteThresh)
    samples = self.samples
    for variant in variants:
      records = [ rec for x in self.vcfReaders for rec in  safe_records(x,variant) if rec] 
      if records:
          ensemble = variant_ensemble(recordSet=[ x for x in records if x], samples=samples, threshold = genoThresh, ignoreMissing = self.ignoreMissing)
          yield records, ensemble.set_concordance()


def safe_records(reader,variant):
    '''
    Returns only records matching the input variant and when there is a ValueError it
    returns [None] this allows it to work in cases where one of the files is empty or
    it doesn't contain some of the variants.
    '''
    try:
        records = [rec for rec in reader.fetch(variant.contig,variant.pos-1,variant.pos) if (len(rec.ALT) == 1 and simple_variant(rec,rec.ALT[0])==variant) ]
    except ValueError:
        records = [None]
    return records


def get_variant_iterator(reader):
    ''' Returns an iterator with sorted variants.
        The variants are trully sorted if the vcf standard is
        followed and there is only one ref for each pos.
    '''
    for rec in reader:
        for alt in sorted(rec.ALT):
            yield simple_variant(rec,alt) 



