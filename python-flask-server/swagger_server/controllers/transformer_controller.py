import connexion
import six

from swagger_server.models.gene_info import GeneInfo  # noqa: E501
from swagger_server.models.transformer_info import TransformerInfo  # noqa: E501
from swagger_server.models.transformer_query import TransformerQuery  # noqa: E501
from swagger_server import util
from swagger_server.models.parameter import Parameter
from swagger_server.models.gene_info import GeneInfoIdentifiers
from swagger_server.models.attribute import Attribute

#REQUIREMENTS
import scipy.stats
from numpy import array, empty
#available at http://software.broadinstitute.org/gsea/downloads.jsp
msigdb_gmt_files=['dat/c2.all.current.0.entrez.gmt', 'dat/c5.all.current.0.entrez.gmt']

#END REQUIREMENTS

transformer_name = 'MSigDB hypergeometric enrichment filter'
valid_controls = ['max p-value', 'max q-value']
default_control_values = {'max p-value': 1e-5, 'max q-value': 0.05}
default_control_types = {'max p-value': 'double', 'max q-value': 'double'}

def get_control(controls, control):
    value = controls[control] if control in controls else default_control_values[control]
    if default_control_types[control] == 'double':
        return float(value)
    elif default_control_types[control] == 'Boolean':
        return bool(value)
    elif default_control_types[control] == 'int':
        return int(value)
    else:
        return value

def entrez_gene_id(gene: GeneInfo):
    """
        Return value of the entrez_gene_id attribute
    """
    if (gene.identifiers is not None and gene.identifiers.entrez is not None):
        if (gene.identifiers.entrez.startswith('NCBIGene:')):
            return gene.identifiers.entrez[9:]
        else:
            return gene.identifiers.entrez
    return None

def correct_pvalues_for_multiple_testing(pvalues, correction_type = "Benjamini-Hochberg"):
    """
    consistent with R - print correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1])
    """
    pvalues = array(pvalues)
    n = int(pvalues.shape[0])
    new_pvalues = empty(n)
    if correction_type == "Bonferroni":
        new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            new_pvalues[i] = (n-rank) * pvalue
    elif correction_type == "Benjamini-Hochberg":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = n - i
            pvalue, index = vals
            new_values.append((n/rank) * pvalue)
        for i in range(0, int(n)-1):
            if new_values[i] < new_values[i+1]:
                new_values[i+1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]
    return new_pvalues

def transform_post(query):  # noqa: E501
    """transform_post

     # noqa: E501

    :param query: Performs transformer query.
    :type query: dict | bytes

    :rtype: List[GeneInfo]
    """
    if connexion.request.is_json:
        query = TransformerQuery.from_dict(connexion.request.get_json())  # noqa: E501

    controls = {control.name:control.value for control in query.controls}
    #Add the originally input genes
    gene_list = query.genes
    genes = dict([(entrez_gene_id(gene_list[i]) if entrez_gene_id(gene_list[i]) != None else gene_list[i].gene_id, i) for i in range(len(gene_list))])

    #Read in the gene sets
    gene_set_y_gene_list_y = {}
    gene_set_y_gene_list_n = {}
    gene_set_n_gene_list_y = {}
    gene_set_n_gene_list_n = {}
    gene_set_k = {}
    gene_set_N = {}
    gene_set_gene_list_gene_ids = {}
    all_gene_set_gene_ids = set()
    for msigdb_gmt_file in msigdb_gmt_files:
        msigdb_gmt_fh = open(msigdb_gmt_file)
        for line in msigdb_gmt_fh:
            cols = line.strip().split('\t')
            if len(cols) < 3:
                continue
            gene_set_id = cols[0]
            gene_ids = cols[2:len(cols)]
            gene_set_gene_list_gene_ids[gene_set_id] = [x for x in gene_ids if x in genes]
            overlap = len(gene_set_gene_list_gene_ids[gene_set_id])
            if overlap == 0:
                continue
            gene_set_y_gene_list_y[gene_set_id] = overlap
            gene_set_N[gene_set_id] = len(gene_ids)

            gene_set_y_gene_list_n[gene_set_id] = gene_set_N[gene_set_id] - gene_set_y_gene_list_y[gene_set_id]
            gene_set_n_gene_list_y[gene_set_id] = len(genes) - gene_set_y_gene_list_y[gene_set_id]
            for x in gene_ids:
                all_gene_set_gene_ids.add(x)
        msigdb_gmt_fh.close()
    M = len(all_gene_set_gene_ids)

    gene_set_pvalues = {}
    gene_set_qvalues = {}
    gene_set_odds_ratios = {}
    all_pvalues = []
    all_gene_set_ids = []
    for gene_set_id in gene_set_y_gene_list_y:
        gene_set_n_gene_list_n[gene_set_id] = M - gene_set_y_gene_list_y[gene_set_id] - gene_set_y_gene_list_n[gene_set_id] - gene_set_n_gene_list_y[gene_set_id]

        table = [[gene_set_y_gene_list_y[gene_set_id], gene_set_y_gene_list_n[gene_set_id]], [gene_set_n_gene_list_y[gene_set_id], gene_set_n_gene_list_n[gene_set_id]]]
        odds_ratio, pvalue = scipy.stats.fisher_exact(table)

        all_pvalues.append(pvalue)
        all_gene_set_ids.append(gene_set_id)

        if pvalue < get_control(controls, 'max p-value'):
            gene_set_pvalues[gene_set_id] = pvalue
            gene_set_odds_ratios[gene_set_id] = odds_ratio

    all_qvalues = correct_pvalues_for_multiple_testing(all_pvalues, correction_type="Benjamini-Hochberg")
    for i, gene_set_id in enumerate(all_gene_set_ids):
        if gene_set_id in gene_set_pvalues and all_qvalues[i] < get_control(controls, 'max q-value'):
            gene_set_qvalues[gene_set_id] = all_qvalues[i]

    final_genes = {}
    final_gene_list = []
    for gene_set_id in sorted(gene_set_qvalues.keys(), key=lambda x: gene_set_qvalues[x]):
        for gene_id in gene_set_gene_list_gene_ids[gene_set_id]:
            gene_entrez_id = "NCBIGene:%s" % gene_id
            if gene_id in genes and gene_entrez_id not in final_genes:
                final_genes[gene_entrez_id] = gene_list[genes[gene_id]]
                final_genes[gene_entrez_id].gene_id = gene_entrez_id
                final_genes[gene_entrez_id].attributes.append(
                      Attribute(
                        name = 'gene set',
                        value = gene_set_id,
                        source = transformer_name
                      ))
                final_genes[gene_entrez_id].attributes.append(
                      Attribute(
                        name = 'p-value',
                        value = gene_set_pvalues[gene_set_id],
                        source = transformer_name
                      ))
                final_genes[gene_entrez_id].attributes.append(
                      Attribute(
                        name = 'q-value',
                        value = gene_set_qvalues[gene_set_id],
                        source = transformer_name
                      ))
                final_genes[gene_entrez_id].attributes.append(
                      Attribute(
                        name = 'odds ratio',
                        value = gene_set_odds_ratios[gene_set_id],
                        source = transformer_name
                      ))
                final_gene_list.append(final_genes[gene_entrez_id])
    return final_gene_list


def transformer_info_get():  # noqa: E501
    """Retrieve transformer info

    Provides information about the transformer. # noqa: E501


    :rtype: TransformerInfo
    """
    return TransformerInfo(
        name = transformer_name,
        function = 'filter',
        #operation = 'enrichment',
        #ui_label = 'HyperGeomEnrich',
        #source_url = 'http://software.broadinstitute.org/gsea/downloads.jsp',
        description = 'Gene-list filter based on hypergeometric enrichment in MSigDB gene sets (http://software.broadinstitute.org/gsea/index.jsp).',
        parameters = [Parameter(x, default_control_types[x], default_control_values[x]) for x in valid_controls],
        required_attributes = ['identifiers.entrez','gene_symbol']
    )

