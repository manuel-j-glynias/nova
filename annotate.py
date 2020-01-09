import json
import time
import requests
import re
import liftover
import xlrd

def get_strand_for_gene(gene,strands):
    strand = '+'
    if strands[gene] == 'Reverse':
        strand = '-'
    return strand

def read_strands(file_name):
    strands = {}
    wb = xlrd.open_workbook(filename=file_name)
    sheet = wb.sheet_by_index(0)
    for row_idx in range(1, sheet.nrows):
        gene = sheet.cell_value(row_idx, 0)
        strand = sheet.cell_value(row_idx, 1)
        strands.update({gene : strand})

    return strands

def add_ckb_gene_info(variant, db):
    gene = variant['gene']
    mycol = db["gene_info"]
    myquery = {'gene': gene}
    mydoc = mycol.find_one(myquery)
    variant['gene_description'] = ''
    variant['gene_category'] = ''
    if mydoc is not None:
        if 'description' in mydoc:
            variant['gene_description'] = mydoc['description']
        if 'category' in mydoc:
            variant['gene_category'] = mydoc['category']

def get_pos_from_pdot(pDot):
    pos = 1
    p = re.compile('\d+\D')
    m = p.search(pDot)
    if m:
        pos = int(m.group()[:-1])
    return pos

def add_ckb_variant_info(variant, db):
    gene = variant['gene']
    pdot = variant['pdot']
    full_name = gene + ' ' + pdot
    mycol = db["variant_info"]
    myquery = {'full_name': full_name}
    mydoc = mycol.find_one(myquery)
    variant['ckb_id'] = ''
    variant['variant_description'] = ''
    variant['protein_effect'] = ''
    variant['variant_type'] = ''
    variant['gDot'] = ''
    variant['cDot'] = ''
    variant['cdot_pos'] = ''
    variant['pdot_pos'] = get_pos_from_pdot(pdot)
    variant['gene_id'] = ''
    variant['is_gain_of_function'] = False
    variant['is_loss_of_function'] = False
    if mydoc is not None:
        if 'ckb_id' in mydoc:
            variant['ckb_id'] = mydoc['ckb_id']
        if 'description' in mydoc:
            variant['variant_description'] = mydoc['description']
        if 'protein_effect' in mydoc:
            variant['protein_effect'] = mydoc['protein_effect']
            variant['is_gain_of_function'] = 'gain of function' in variant['protein_effect']
            variant['is_loss_of_function'] = 'loss of function' in variant['protein_effect']
        if 'variant_type' in mydoc:
            variant['variant_type'] = mydoc['variant_type']
        if 'gDot' in mydoc:
            variant['gDot'] = mydoc['gDot']
        if 'cDot' in mydoc:
            variant['cDot'] = mydoc['cDot']
        if 'cdot_pos' in mydoc:
            variant['cdot_pos'] = mydoc['cdot_pos']
        if 'pdot_pos' in mydoc:
            variant['pdot_pos'] = mydoc['pdot_pos']
        if 'gene_id' in mydoc:
            variant['gene_id'] = mydoc['gene_id']


def get_num_cv_nonbenign_between(gene, begin, end, db):
    mycol = db["clinvar"]
    myquery = {'gene': gene, 'pdot_pos': {'$gte': begin, '$lte': end},
               'is_majority_vote_not_benign': True}
    count = mycol.find(myquery).count()
    return count


def get_num_cv_pathogenic_between(gene, begin, end, db):
    mycol = db["clinvar"]
    myquery = {'gene': gene, 'pdot_pos': {'$gte': begin, '$lte': end},
               'is_majority_vote_pathogenic': True}
    count = mycol.find(myquery).count()
    return count


def has_cv_pathogenic_or_LOF_variants_on_both_sides(gene, pos, window, db):
    b = False
    n = get_num_cv_pathogenic_between(gene, pos - window, pos, db)
    if n == 0:
        n = get_num_ckb_lof_variants_between(gene, pos - window, pos, db)
    m = get_num_cv_pathogenic_between(gene, pos, pos + window, db)
    if m == 0:
        m = get_num_ckb_lof_variants_between(gene, pos, pos + window, db)
    if n > 0 and m > 0:
        b = True
    return b, n, m


def get_num_ckb_lof_variants_between(gene, begin, end, db):
    mycol = db["variant_info"]
    myquery = {'gene': gene, 'pdot_pos': {'$gte': begin, '$lte': end},
               '$or': [{'protein_effect': 'loss of function'}, {'protein_effect': 'loss of function - predicted'}]}
    count = mycol.find(myquery).count()
    return count


def get_num_ckb_gof_variants_between(gene, begin, end, db):
    mycol = db["variant_info"]
    myquery = {'gene': gene, 'pdot_pos': {'$gte': begin, '$lte': end},
               '$or': [{'protein_effect': 'gain of function'}, {'protein_effect': 'gain of function - predicted'}]}
    count = mycol.find(myquery).count()
    return count


def add_hot_spot_info(variant,db):
    gene = variant['gene']
    mycol = db["hotspots"]
    myquery = {'gene': gene}
    hotspots = []
    mydocs = mycol.find(myquery).sort("start")
    if mydocs is not None:
        for doc in mydocs:
            region = {}
            region['begin'] = doc['start']
            region['end'] = doc['end']
            region['residue'] = doc['residue']
            # region['ckb_count'] = get_num_ckb_gof_variants_between(gene, region['begin'], region['end'])
            # region['cv_count'] = get_num_cv_nonbenign_between(gene, region['begin'], region['end'])
            hotspots.append(region)
    variant['hotspots'] = hotspots

def add_clinvar(variant, db):
    gene = variant['gene']
    pdot = variant['pdot']
    mycol = db["clinvar"]
    myquery = {'gene': gene, 'pDot': pdot}
    mydoc = mycol.find_one(myquery)
    if (mydoc != None):
        variant['in_clinvar'] = True
        variant['clinvar_significance'] = mydoc['significance']
        variant['clinvar_explain'] = mydoc['explain']
        variant['variant_id'] = mydoc['variant_id']
        variant['is_clinvar_not_benign'] = mydoc['is_majority_vote_not_benign']
        variant['is_clinvar_benign'] = mydoc['is_majority_vote_benign']
        variant['is_clinvar_pathogenic'] = mydoc['is_majority_vote_pathogenic']
    else:
        variant['in_clinvar'] = False


def is_near_LOF_mutation(variant, db):
    window = 10
    gene = variant['gene']
    pos = variant['pdot_pos']
    if (gene != None and pos != None and pos != ''):
        b, upstream, downstream = has_cv_pathogenic_or_LOF_variants_on_both_sides(gene, int(pos), window,db)
        variant['is_near_LOF_mutation'] = b
        variant['num_upstream_LOF_mutations'] = upstream
        variant['num_downstream_LOF_mutations'] = downstream
    else:
        variant['is_near_LOF_mutation'] = False


def is_better(appris, best_appris):
    appris_val = int(appris[1:])
    best_val = int(best_appris[1:])
    return appris_val < best_val

def annotate_one(variant,db):
    protein_altering_variants = ['transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant',
                                 'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost',
                                 'transcript_amplification',
                                 'inframe_insertion', 'inframe_deletion', 'missense_variant',
                                 'protein_altering_variant']
    truncating_variants = ['transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained',
                           'frameshift_variant', 'start_lost']
    response_dict = {'status': 'unspecified parameters'}

    server = "https://rest.ensembl.org/vep/human/region/"
    ext = str(variant['chrom']) + ':' + str(variant['pos_hg38']) + '/' + variant['alt_hg38'] + '?appris=1&'
    # print(server + ext)
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    # if variant['chrom']=='17':
    #     print('17')
    if not r.ok:
        #         r.raise_for_status()
        print(r.text)
        response_dict['status'] = 'server error'
    else:
        decoded = r.json()
        best_appris = 'R9'
        if 'transcript_consequences' in decoded[0]:
            for consequence in decoded[0]['transcript_consequences']:
                if 'protein_start' in consequence and 'appris' in consequence and consequence['appris'].startswith(\
                        'P'):
                    if is_better(consequence['appris'], best_appris):
                        best_appris = consequence['appris']
                        variant['gene'] = consequence['gene_symbol']
                        add_ckb_gene_info(variant,db)
                        # add_uniprot_info(variant)
                        mutation_type = consequence['consequence_terms'][0]
                        if 'cds_start' in consequence:
                            allele_string = str(decoded[0]['allele_string'])
                            if '-' in allele_string:
                                allele_string = 'del'
                            else:
                                allele_string = allele_string.replace('/','>')
                            variant['cdot'] = str(consequence['cds_start']) + allele_string
                        if 'protein_start' in consequence:
                            protein_start = consequence['protein_start']
                            variant['protein_start'] = protein_start
                            # add_region_hit(variant)
                            # add_hot_spot_info(variant)

                            if 'amino_acids' in consequence:
                                amino_acids = consequence['amino_acids']
                                aa = amino_acids.split('/')
                                variant['pdot'] = aa[0] + str(protein_start) + aa[1]
                            else:
                                variant['pdot'] = str(protein_start) + 'del'
                            add_ckb_variant_info(variant,db)
                            variant['full_name'] = variant['gene'] + ' ' + variant['pdot']
                        else:
                            print('crap')
                        variant['mutation_type'] = mutation_type
                        variant['is_protein_altering'] = mutation_type in protein_altering_variants
                        variant['is_truncating_variants'] = mutation_type in truncating_variants
                        if 'polyphen_prediction' in consequence:
                            variant['polyphen_prediction'] = consequence['polyphen_prediction']
                        if 'sift_prediction' in consequence:
                            variant['sift_prediction'] = consequence['sift_prediction']
                        variant['predicted_deleterious'] = False
                        if 'polyphen_prediction' in variant and 'sift_prediction' in variant:
                            if 'damaging' in variant['polyphen_prediction'] and 'deleterious' in variant[
                                'sift_prediction']:
                                variant['predicted_deleterious'] = True

                        response_dict['status'] = 'success'
                        response_dict['variant'] = variant
    return response_dict

def annotate_one_by_pick(variant,db):
    protein_altering_variants = ['transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant',
                                 'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost',
                                 'transcript_amplification',
                                 'inframe_insertion', 'inframe_deletion', 'missense_variant',
                                 'protein_altering_variant']
    truncating_variants = ['transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained',
                           'frameshift_variant', 'start_lost']
    response_dict = {'status': 'unspecified parameters'}

    server = "https://rest.ensembl.org/vep/human/region/"
    ext = str(variant['chrom']) + ':' + str(variant['pos_hg38']) + '/' + variant['alt_hg38'] + '?pick=1'
    # print(server + ext)
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    # if variant['chrom']=='17':
    #     print('17')
    if not r.ok:
        #         r.raise_for_status()
        print(r.text)
        response_dict['status'] = 'server error'
    else:
        decoded = r.json()
        if 'transcript_consequences' in decoded[0]:
            for consequence in decoded[0]['transcript_consequences']:
                if 'protein_start' in consequence :
                        variant['gene'] = consequence['gene_symbol']
                        add_ckb_gene_info(variant,db)
                        # add_uniprot_info(variant)
                        mutation_type = consequence['consequence_terms'][0]
                        if 'cds_start' in consequence:
                            allele_string = str(decoded[0]['allele_string'])
                            if '-' in allele_string:
                                allele_string = 'del'
                            else:
                                allele_string = allele_string.replace('/','>')
                            variant['cdot'] = str(consequence['cds_start']) + allele_string
                        if 'protein_start' in consequence:
                            protein_start = consequence['protein_start']
                            variant['protein_start'] = protein_start
                            # add_region_hit(variant)
                            # add_hot_spot_info(variant)

                            if 'amino_acids' in consequence:
                                amino_acids = consequence['amino_acids']
                                aa = amino_acids.split('/')
                                variant['pdot'] = aa[0] + str(protein_start) + aa[1]
                            else:
                                variant['pdot'] = str(protein_start) + 'del'
                            add_ckb_variant_info(variant,db)
                            variant['full_name'] = variant['gene'] + ' ' + variant['pdot']
                        else:
                            print('crap')
                        variant['mutation_type'] = mutation_type
                        variant['is_protein_altering'] = mutation_type in protein_altering_variants
                        variant['is_truncating_variants'] = mutation_type in truncating_variants
                        if 'polyphen_prediction' in consequence:
                            variant['polyphen_prediction'] = consequence['polyphen_prediction']
                        if 'sift_prediction' in consequence:
                            variant['sift_prediction'] = consequence['sift_prediction']
                        variant['predicted_deleterious'] = False
                        if 'polyphen_prediction' in variant and 'sift_prediction' in variant:
                            if 'damaging' in variant['polyphen_prediction'] and 'deleterious' in variant[
                                'sift_prediction']:
                                variant['predicted_deleterious'] = True

                        response_dict['status'] = 'success'
                        response_dict['variant'] = variant
    return response_dict


def call_annotate(variant,db):
    variant['is_protein_altering'] = False
    variant['predicted_deleterious'] = False
    response_dict = annotate_one(variant, db)
    if response_dict['status'] == 'success':
        variant = response_dict['variant']
        add_clinvar(variant,db)
        add_hot_spot_info(variant,db)
        is_near_LOF_mutation(variant,db)

    else:
        print('annotate_one failed with response status:',response_dict['status'])
        response_dict = annotate_one_by_pick(variant, db)
        if response_dict['status'] == 'success':
            variant = response_dict['variant']
            add_clinvar(variant,db)
            add_hot_spot_info(variant,db)
            is_near_LOF_mutation(variant,db)
        else:
            print('annotate_one_by_pick failed with response status:', response_dict['status'])



def get_annotated_snv(variant, db):
    ann_var = None
    mycol = db["ann_var"]
    key = get_key_for_variant(variant)
    myquery = {'key': key}
    mydoc = mycol.find_one(myquery)
    if mydoc is not None:
        json_string = mydoc['value']
        ann_var = json.loads(json_string)
    return ann_var


def get_key_for_variant(variant):
    key = variant['chrom'] + '_' + str(variant['pos_hg19']) + '_' + variant['alt']
    return key

def save_annotated_snv(variant, db):
    key = get_key_for_variant(variant)
    json_string = json.dumps(variant)
    mycol = db["ann_var"]
    kv_dict = {'key': key, 'value':json_string}
    mycol.insert_one(kv_dict)

def handle_liftover(strands, variant):
    variant['pos_hg38'] = liftover.convert_coordinated(variant['chrom'], variant['pos_hg19'])
    variant['strand'] = get_strand_for_gene(variant['HGNC_Symbol'], strands)
    variant['alt_hg38'] = variant['alt']
    if len(variant['ref']) > len(variant['alt_hg38']):
        #     is deletion
        if variant['strand'] == '+':
            if len(variant['ref']) - len(variant['alt_hg38']) == 1:
                variant['pos_hg38'] = variant['pos_hg38'] + 1
                variant['alt_hg38'] = '-'
            else:
                del_length = len(variant['ref']) - len(variant['alt_hg38'])
                variant['pos_hg38'] = str(variant['pos_hg38'] + 1) + '..' + str(variant['pos_hg38'] + del_length)
                variant['alt_hg38'] = '-'
        else:
            if len(variant['ref']) - len(variant['alt_hg38']) == 1:
                variant['pos_hg38'] = variant['pos_hg38'] - 1
                variant['alt_hg38'] = '-'
            else:
                del_length = len(variant['ref']) - len(variant['alt_hg38'])
                pos = variant['pos_hg38'] - 1
                # variant['pos_hg38'] = str(pos) + '..' + str(pos - del_length)
                variant['pos_hg38'] = str(pos - del_length) + '..' + str(pos)
                variant['alt_hg38'] = '-'

    elif len(variant['ref']) < len(variant['alt_hg38']):
        #       is insertion
        # variant['pos_hg38'] = str(variant['pos_hg38']) + '..' + str(variant['pos_hg38'])
        if len(variant['ref'])==1:
            if variant['strand'] == '+':
                variant['pos_hg38'] = str(variant['pos_hg38']+1)
            else:
                variant['pos_hg38'] = str(variant['pos_hg38'] - 1)


def annotate_snv(strands,variant,db):
    annotated_snv = get_annotated_snv(variant,db)
    if annotated_snv is None:
        handle_liftover(strands, variant)
        call_annotate(variant, db)
        save_annotated_snv(variant,db)
        annotated_snv = variant
    else:
        annotated_snv['VAF'] = variant['VAF']
        annotated_snv['Quality'] = variant['Quality']
    return annotated_snv

def annotate_cnv(variant,db):
    add_ckb_gene_info(variant, db)
    return variant