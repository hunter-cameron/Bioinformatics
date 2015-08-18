from html.parser import HTMLParser

class Parser(HTMLParser):

    def __init__(self):
        super(Parser, self).__init__()
        self.data_dict = {}
    
        # stores the current data key
        self.cur_key = None

        # stores what to do with the data
        self.slurp_mode = None


    def handle_starttag(self, tag, attr):
        if tag == "th":
            # headers are the keys
            self.slurp_mode = "key"
        elif tag == "td":
            self.slurp_mode = "data"
        
    def handle_data(self, data):
        data = data.strip()
        if data:
            if self.slurp_mode == "key":
                self.cur_key = data.strip()
                print(self.cur_key)
            elif self.slurp_mode == "data":
                self.data_dict[self.cur_key] = self.data_dict.get(self.cur_key, "") + data.strip()







data = """
<a name='overview' href='#'><h2>Overview</h2> </a>
<table class='img'  border='1' >
<tr class='img' >
  <th class='subhead' align='right'>Study Name (Proposal Name)</th>
  <td class='img'   align='left'>Rhizosphere Grand Challenge Isolate Sequencing</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Organism Name</th>
  <td class='img'   align='left'>Arthrobacter sp. 131MFCol6.1</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Taxon ID</th>
  <td class='img'   align='left'>2517572122</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>IMG Submission ID</th>
  <td class='img'   align='left'><a href='https://img.jgi.doe.gov/cgi-bin/submit/main.cgi?section=ERSubmission&page=displaySubmission&submission_id=10932'  onclick="">10932</a></td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>NCBI Taxon ID</th>
  <td class='img'   align='left'><a href='http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1157944'  onclick="">1157944</a></td>
</tr>
<tr class='img'>
<th class='subhead'>
GOLD ID in IMG Database</th>
</td>
<td class='img'>
<a href='https://gold.jgi-psf.org/study?id=Gs0030678'  onclick="">Study ID: Gs0030678</a>&nbsp;&nbsp;<a href='https://gold.jgi-psf.org/projects?id=Gp0018766'  onclick="">Project ID: Gp0018766</a></td>
</tr>
<tr class='img'>
<th class='subhead'>
GOLD Analysis Project Id</th>
</td>
<td class='img'>
<a href='https://gold.jgi-psf.org/analysis_projects?id=Ga0000025'  onclick="">Ga0000025</a></td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>GOLD Analysis Project Type</th>
  <td class='img'   align='left'>Genome Analysis</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Submission Type</th>
  <td class='img'   align='left'>Primary</td>
</tr>
<tr class='img'>
<th class='subhead'>External Links</th>
<td class='img'>
<a href="http://genome.jgi.doe.gov/lookup?keyName=jgiProjectId&keyValue=1000073" onClick="_gaq.push(['_trackEvent', 'Download Data', 'JGI Portal', '2517572122']);" > 
               <img style='border:none; vertical-align:text-top;' 
                 src='https://img.jgi.doe.gov/w/images/genomeProjects_icon.gif' />
             JGI Portal </a>&nbsp; </td>
</tr>
<tr class='img' >
<th class='subhead' align='right'>Lineage</th>
<td class='img' ><a href='main.cgi?section=TaxonList&page=lineageMicrobes&domain=Bacteria'  onclick="">Bacteria</a>; <a href='main.cgi?section=TaxonList&page=lineageMicrobes&phylum=Actinobacteria'  onclick="">Actinobacteria</a>; <a href='main.cgi?section=TaxonList&page=lineageMicrobes&ir_class=Actinobacteria'  onclick="">Actinobacteria</a>; <a href='main.cgi?section=TaxonList&page=lineageMicrobes&ir_order=Actinomycetales'  onclick="">Actinomycetales</a>; <a href='main.cgi?section=TaxonList&page=lineageMicrobes&family=Micrococcaceae'  onclick="">Micrococcaceae</a>; <a href='main.cgi?section=TaxonList&page=lineageMicrobes&genus=Arthrobacter'  onclick="">Arthrobacter</a>; <a href='main.cgi?section=TaxonList&page=lineageMicrobes&species=Arthrobacter+sp.+131MFCol6.1'  onclick="">Arthrobacter sp. 131MFCol6.1</a></td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Sequencing Status</th>
  <td class='img'   align='left'>Permanent Draft</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Sequencing Center</th>
  <td class='img'   align='left'>DOE Joint Genome Institute (JGI)</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>IMG Release</th>
  <td class='img'   align='left'>&nbsp; </td>
</tr>
<tr class='img' >
<th class='subhead'>Comment</th>
<td class='img' >
&nbsp; 
</td>
</tr>
<tr class='img' >
<th class='subhead'> Release Date </th>
<td class='img' >
2013-05-14</td></tr>
<tr class='img' >
<th class='subhead'> Add Date </th>
<td class='img' >
2012-10-31</td></tr>
<tr class='img' >
<th class='subhead'> Modified Date </th>
<td class='img' >
2014-08-05</td></tr>
<tr class='img' >
<th class='subhead'> Distance Matrix Calc. Date </th>
<td class='img' >
2014-07-02</td></tr>
<tr class='img' >
<th class='subhead'> High Quality </th>
<td class='img' >
Yes</td></tr>
<tr class='img' >
<th class='subhead'>IMG Product Flag </th>
<td class='img' >
Yes</td></tr>
<tr class='img' >
  <th class='subhead' align='right'>Is Public</th>
  <td class='img'   align='left'>Yes</td>
</tr>
<tr class='highlight'>
<th class='subhead'>Project Information</th>  <th class='subhead'> &nbsp; </th></tr>
<tr class='img' >
  <th class='subhead' align='right'>Bioproject Accession</th>
  <td class='img'   align='left'>PRJNA175167</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Biosample Accession</th>
  <td class='img'   align='left'>SAMN02441022</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Culture Type</th>
  <td class='img'   align='left'>Isolate</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>GOLD Sequencing Strategy</th>
  <td class='img'   align='left'>Whole Genome Sequencing</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Gram Staining</th>
  <td class='img'   align='left'>Gram+</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Seq Status</th>
  <td class='img'   align='left'>Complete</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Sequencing Method</th>
  <td class='img'   align='left'>Illumina, Illumina HiSeq 2000, Illumina HiSeq 2500</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Type Strain</th>
  <td class='img'   align='left'>Unknown</td>
</tr>
</tr>
<tr class='highlight'>
<th class='subhead'>Phenotypes/Metabolism from Pathway Assertion</th> <th class='subhead'> &nbsp; </th> 
<tr class='img' >
  <th class='subhead' align='right'>Metabolism</th>
  <td class='img'   align='left'><a href='main.cgi?section=TaxonDetail&page=taxonPhenoRuleDetail&taxon_oid=2517572122&rule_id=7'  onclick="">Auxotroph (L-lysine auxotroph)</a> (IMG_PIPELINE; 2013-12-02)</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Metabolism</th>
  <td class='img'   align='left'><a href='main.cgi?section=TaxonDetail&page=taxonPhenoRuleDetail&taxon_oid=2517572122&rule_id=12'  onclick="">Prototrophic (L-aspartate prototroph)</a> (IMG_PIPELINE; 2013-12-02)</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Metabolism</th>
  <td class='img'   align='left'><a href='main.cgi?section=TaxonDetail&page=taxonPhenoRuleDetail&taxon_oid=2517572122&rule_id=14'  onclick="">Prototrophic (L-glutamate prototroph)</a> (IMG_PIPELINE; 2013-12-02)</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Metabolism</th>
  <td class='img'   align='left'><a href='main.cgi?section=TaxonDetail&page=taxonPhenoRuleDetail&taxon_oid=2517572122&rule_id=17'  onclick="">Auxotroph (L-phenylalanine auxotroph)</a> (IMG_PIPELINE; 2013-12-02)</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Metabolism</th>
  <td class='img'   align='left'><a href='main.cgi?section=TaxonDetail&page=taxonPhenoRuleDetail&taxon_oid=2517572122&rule_id=18'  onclick="">Auxotroph (L-tyrosine auxotroph)</a> (IMG_PIPELINE; 2013-12-02)</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Metabolism</th>
  <td class='img'   align='left'><a href='main.cgi?section=TaxonDetail&page=taxonPhenoRuleDetail&taxon_oid=2517572122&rule_id=21'  onclick="">Auxotroph (L-tryptophan auxotroph)</a> (IMG_PIPELINE; 2013-12-02)</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Metabolism</th>
  <td class='img'   align='left'><a href='main.cgi?section=TaxonDetail&page=taxonPhenoRuleDetail&taxon_oid=2517572122&rule_id=22'  onclick="">Auxotroph (L-histidine auxotroph)</a> (IMG_PIPELINE; 2013-12-02)</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Metabolism</th>
  <td class='img'   align='left'><a href='main.cgi?section=TaxonDetail&page=taxonPhenoRuleDetail&taxon_oid=2517572122&rule_id=23'  onclick="">Prototrophic (Glycine prototroph)</a> (IMG_PIPELINE; 2013-12-02)</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Metabolism</th>
  <td class='img'   align='left'><a href='main.cgi?section=TaxonDetail&page=taxonPhenoRuleDetail&taxon_oid=2517572122&rule_id=26'  onclick="">Auxotroph (L-arginine auxotroph)</a> (IMG_PIPELINE; 2013-12-02)</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Metabolism</th>
  <td class='img'   align='left'><a href='main.cgi?section=TaxonDetail&page=taxonPhenoRuleDetail&taxon_oid=2517572122&rule_id=34'  onclick="">Auxotroph (L-isoleucine auxotroph)</a> (IMG_PIPELINE; 2013-12-02)</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Metabolism</th>
  <td class='img'   align='left'><a href='main.cgi?section=TaxonDetail&page=taxonPhenoRuleDetail&taxon_oid=2517572122&rule_id=36'  onclick="">Auxotroph (L-leucine auxotroph)</a> (IMG_PIPELINE; 2013-12-02)</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Metabolism</th>
  <td class='img'   align='left'><a href='main.cgi?section=TaxonDetail&page=taxonPhenoRuleDetail&taxon_oid=2517572122&rule_id=44'  onclick="">Auxotroph (L-valine auxotroph)</a> (IMG_PIPELINE; 2013-12-02)</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Metabolism</th>
  <td class='img'   align='left'><a href='main.cgi?section=TaxonDetail&page=taxonPhenoRuleDetail&taxon_oid=2517572122&rule_id=46'  onclick=""> (Non-selenocysteine synthesizer)</a> (IMG_PIPELINE; 2013-12-02)</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Metabolism</th>
  <td class='img'   align='left'><a href='main.cgi?section=TaxonDetail&page=taxonPhenoRuleDetail&taxon_oid=2517572122&rule_id=48'  onclick=""> (Non-biotin synthesizer)</a> (IMG_PIPELINE; 2013-12-02)</td>
</tr>
<tr class='img' >
  <th class='subhead' align='right'>Metabolism</th>
  <td class='img'   align='left'><a href='main.cgi?section=TaxonDetail&page=taxonPhenoRuleDetail&taxon_oid=2517572122&rule_id=50'  onclick="">Auxotroph (Incomplete Coenyzme A biosynthesis)</a> (IMG_PIPELINE; 2013-12-02)</td>
</tr>
</tr>
</table>
"""


parser = Parser()
parser.feed(data)


for k, v in parser.data_dict.items():
    print((k, v))
